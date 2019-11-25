classdef KFESolver
	% This is a solver for the Kolmogorov-Forward Equation, providing
	% both a direct solver and an iterative solver.

	properties (SetAccess=protected)
		p;
		income;
		grdKFE;
		options;
        n_states;
	end

	methods
		function obj = KFESolver(p, income, grdKFE)
			% Parameters
			% ----------

			% p : 	An object with the following required attributes:
			%		nb_KFE (number of points on liquid asset grid)
			%		na_KFE (number of points on illiquid asset grid)
			%		nz (number of states in extra dimension of heterogeneity)
			%		deathrate (rate of death)
			%		ResetIncomeUponDeath (whether or not income is redrawn upon death, untested...)
			%
			%		and (optionally) any of the following attributes:
			%		deltaKFE (default 1e5, step size)
			%		KFE_tol (default 1e-8)
            %       KFE_maxiter (default 1e-4)
			%		iterateKFE (default true, use iterative rather than direct method)
			%		intermediateCheck (default true, throws error if convergence appears unlikely)

			% income : An object with the following required attributes:
			%		ny (number of income grid points)
			%		ytrans (the income transition matrix, row sums should be 0)
			%		ydist (stationary pmf of the income process)
				
			% grdKFE : A Grid object
			
			obj.p = p;
			obj.income = income;
			obj.grdKFE = grdKFE;
			obj.n_states = p.nb_KFE * p.na_KFE * p.nz * income.ny;

			required_parameter_vars = {'nb_KFE', 'na_KFE', 'nz',...
					'deathrate', 'ResetIncomeUponDeath'};
			aux.check_for_required_properties(p, required_parameter_vars);

			required_income_vars = {'ny', 'ytrans', 'ydist'};
			aux.check_for_required_properties(income, required_income_vars);

			obj.options = obj.parse_input(aux.to_structure(p));
		end

		function g = solve(obj, A, g0)
			% Parameters
			% ----------

			% A : square, sparse transition matrix
			
			% g0 : optional, the initial distribution for
			%	the iterative method

			% Returns
			% -------
			% g : the stationary distribution

			size_A = size(A);
			assert(size_A(1) == size_A(2), "KFESolver requires a square transition matrix")
			assert(numel(size_A) == 2, "KFESolver requires a square transition matrix")
			assert(issparse(A), "KFESolver requires a sparse transition matrix")

			if obj.options.iterateKFE
				if ~exist('g0')
					% g0 wasn't passed
				    g0 = obj.guess_initial_distribution();
				else
					assert(numel(g) == size_A(1),...
						"Initial distribution and transition matrix have inconsistent size")
				end

			    g = obj.solve_iterative(A, g0);
			else
				g = obj.solve_direct(A);
			end
		end

		function set_option(obj, keyword, value)
			if isprop(keyword)
				obj.options.(keyword) = value;
			end
		end
	end

	methods (Access=private)
		function options = parse_input(obj, p)
			parser = inputParser;
			parser.KeepUnmatched = true;
			addParameter(parser, 'deltaKFE', 1e5);
			addParameter(parser, 'KFE_tol', 1e-8);
			addParameter(parser, 'KFE_maxiter', 1e3);
			addParameter(parser, 'iterateKFE', true);
			addParameter(parser, 'intermediateCheck', true);
			parse(parser, p);
			options = parser.Results;
		end

		function g0 = guess_initial_distribution(obj)
			g0 = ones(obj.p.nb_KFE, obj.p.na_KFE, obj.p.nz, obj.income.ny);
		    g0 = g0 .* permute(repmat(obj.income.ydist,...
		    	[1 obj.p.nb_KFE obj.p.na_KFE obj.p.nz]),[2 3 4 1]);
		    if obj.p.OneAsset
		        g0(:,obj.grdKFE.a.vec>0,:,:) = 0;
		    end
		    g0 = g0 / sum(g0(:));
		    g0 = g0 ./ obj.grdKFE.trapezoidal.matrix;
		end

		function g = solve_direct(obj, A)
			% solves the eigenvalue problem directly

			inctrans = kron(obj.income.ytrans, speye(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz));
	        Ap_extended = sparse([(A+inctrans)'; ones(1, obj.n_states)]);
	        RHS = sparse([zeros(obj.n_states, 1); 1]);

			g = Ap_extended \ RHS;
	        g = g ./ obj.grdKFE.trapezoidal.matrix(:);
	        g = reshape(full(g), obj.p.nb_KFE, obj.p.na_KFE, obj.p.nz, obj.income.ny);
		end

		function g = solve_iterative(obj, A, g0)
			% solves using an iterative method

			% compute the LHS of the KFE
			KFE_LHS = obj.KFE_matrix_divisor(A);

		    % transition matrix with diagonal killed
			ytrans0  = obj.income.ytrans - diag(diag(obj.income.ytrans)); 
			ytrans0p = ytrans0';

			states_per_income = obj.p.nb_KFE * obj.p.na_KFE * obj.p.nz;
			iter = 0;
			dst = 1e5;
			fprintf('    --- Iterating over KFE ---\n')
            g = g0(:);
			while (iter <= obj.options.KFE_maxiter) && (dst > obj.options.KFE_tol)
				iter = iter + 1;

			    gg_tilde = obj.grdKFE.trapezoidal.diagm * g(:);
			    g1 = zeros(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz, obj.income.ny);
			    for iy = 1:obj.income.ny    
			    	gk_sum = sum(repmat(ytrans0p(iy,:), states_per_income, 1)...
			            .* reshape(gg_tilde, states_per_income, obj.income.ny),2);
			    	death_inflows = obj.compute_death_inflows(gg_tilde, iy);
		            g1(:,iy) = KFE_LHS{iy}*(gg_tilde(1+(iy-1)*states_per_income:iy*states_per_income)...
		                		+ obj.options.deltaKFE*gk_sum + obj.options.deltaKFE*death_inflows);
		        end

			    g1 = g1(:) ./ sum(g1(:));
			    g1 = obj.grdKFE.trapezoidal.diagm \ g1;

		        dst = max(abs(g1(:) - g(:)));

		        if obj.options.intermediateCheck
		        	check_if_not_converging(dst, iter);
		        end
		        
			    if (iter==1) || (mod(iter, 100) == 0)
			        fprintf('\tKFE iteration  = %i, distance = %e\n', iter, dst);
			    end
			    g = g1;
			end
			obj.check_if_converged(dst, iter);
			g = reshape(g, obj.p.nb_KFE, obj.p.na_KFE, obj.p.nz, obj.income.ny);
		end

		function LHS = KFE_matrix_divisor(obj, A)
			% Returns
			% -------
			% LHS : a cell array of operators B_k s.t. B_k * RHS_k
			%	returns the k-th income section of the equilibrium distribution,
			%	which is LHS_k \ RHS_k

			states_per_income = obj.p.nb_KFE * obj.p.na_KFE * obj.p.nz;
			LHS = cell(1, obj.income.ny);
			for k = 1:obj.income.ny
				i1 = 1 + (k-1) * states_per_income;
				i2 = k * states_per_income;

				LHS{k} = (speye(states_per_income) - obj.options.deltaKFE * A(i1:i2, i1:i2)'...
			   		- obj.options.deltaKFE * (obj.income.ytrans(k,k) - obj.p.deathrate) * speye(states_per_income));
				LHS{k} = inverse(LHS{k});
			end
		end

		function death_inflows  = compute_death_inflows(obj, gg_tilde, iy)
	    	if (obj.p.Bequests == 1) && (obj.p.ResetIncomeUponDeath == 1)
                death_inflows = obj.p.deathrate * obj.income.ydist(iy) * sum(reshape(gg_tilde, [], obj.income.ny), 2);
            elseif (obj.p.Bequests == 1) && (obj.p.ResetIncomeUponDeath == 0)
                death_inflows = obj.p.deathrate * gg_tilde(1+(iy-1)*(obj.nb_KFE*obj.na_KFE*obj.nz):iy*(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz));
            elseif (obj.p.Bequests == 0) && (obj.p.ResetIncomeUponDeath == 1)
                death_inflows = sparse(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz,1);
                death_inflows(1:obj.p.nb_KFE*obj.p.na_KFE:end) = obj.p.deathrate * obj.income.ydist(iy) * (1/obj.p.nz);
            elseif (obj.p.Bequests == 0) && (obj.p.ResetIncomeUponDeath == 0)
                death_inflows = sparse(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz,1);
                death_inflows(obj.grdKFE.loc0b0a:obj.p.nb_KFE*obj.p.na_KFE:end) = obj.p.deathrate * obj.income.ydist(iy) * (1/obj.p.nz);
	    	end
	    end

	    function check_if_converged(obj, dst, iter)
	    	if dst < obj.options.KFE_tol
			    fprintf('\tKFE converged after %i iterations\n', iter);
			elseif dst >= obj.options.KFE_tol
				error('KFE did not converge')
			end
		end
	end
end

function check_if_not_converging(dst, iter)
	if (dst>10000) && (iter>2000)
    	msgID = 'KFE:NotConverging';
	    msg = 'KFE:NotConverging';
	    KFEException = MException(msgID,msg);
	    throw(KFEException)
    end
end