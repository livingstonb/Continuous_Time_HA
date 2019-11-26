classdef KFESolver
	% This is a solver for the Kolmogorov-Forward Equation, providing
	% both a direct solver and an iterative solver.
	%
	% See the documentation with the 'help KFESolver' command.

	properties (SetAccess=protected)
		%	An object with the following required attributes:
		%
		%		nb_KFE >= 1
		%		- Number of points on the liquid asset grid.
		%
		%		na_KFE >= 1
		%		- Number of points on the illiquid asset grid.
		%
		%		nz >= 1
		%		- Number of states in the extra dimension of
		%		  heterogeneity.
		%
		%		deathrate >= 0
		%		- Poisson rate of death.
		%
		%		ResetIncomeUponDeath, true/false
		%		- Boolean indicator of whether or not income should be
		%		  redrawn from the stationary distribution upon death.
		p;

		%	An object with the following required attributes:
		%
		%		ny >= 1
		%		- number of income grid points
		%
		%		ytrans
		%		- the square income transition matrix, row sums should be 0
		%
		%		ydist
		%		- stationary pmf of the income process, vector, must sum to 1
		income;

		%	A Grid object. See the Grid documentation for details.
		grdKFE;

		%	A KFEOptions object used internally by KFESolver.
		%	The options can be set manually using the set_option()
		%	method, or they can be passed as a KFEOptions object
		%	to the constructor.
		%
		%		delta > 0
		%		- Step size for the iterative method. Default = 1e5.
		%
		%		maxiters > 1
		%		- Max number of iterations for the iterative method.
		%		  Default = 1e4.
		%
		%		tol > 0
		%		- Convergence tolerance for the iterative method.
		%		  Default = 1e-8.
		%
		%		iterate
		%		- Boolean, the solver uses the iterative method if
		%		  iterateKFE is true. Default = true.
		%
		%		intermediate_check
		%		- Boolean, if intermediate_check evaluates to true,
		%		  the solver will check the norm of the error. If
		%		  it is large, an exception is thrown.
		%		  Default = True.
		options;

		% Total number of states.
        n_states;
	end

	methods
		function obj = KFESolver(p, income, grdKFE, options)
			% Class constructor. If p is not a Params object
			% or income is not an Income object, see the help
			% documentation for the requirements these variables
			% must satisfy.
			%
			% The optional argument 'options' must be a KFEOptions
			% object. If this argument is not passed, the default
			% options will be used. See KFEOptions for more details.
			% Note that if 'p' has a KFEOptions object as an
			% attribute named 'kfe_options', that object will
			% will be used to set the options. Passing options
			% to the last argument of this constructor will
			% override this, however.
			
			obj.p = p;
			obj.income = income;
			obj.grdKFE = grdKFE;
			obj.n_states = p.nb_KFE * p.na_KFE * p.nz * income.ny;

			required_parameter_vars = {'nb_KFE', 'na_KFE', 'nz',...
					'deathrate', 'ResetIncomeUponDeath'};
			aux.check_for_required_properties(p, required_parameter_vars);

			required_income_vars = {'ny', 'ytrans', 'ydist'};
			aux.check_for_required_properties(income, required_income_vars);

			options_passed_in_p = false;
			if isprop(p, 'kfe_options')
				if isa(p.kfe_options, 'KFEOptions')
					options_passed_in_p = true;
				end
			end

			if exist('options')
				assert(isa(options, KFEOptions),...
					"options argument must be a KFEOptions object");
				obj.options = options;
			elseif options_passed_in_p
				obj.options = p.kfe_options;
			else
				obj.options = solver.KFEOptions();
			end
		end

		function g = solve(obj, A, g0)
			% Parameters
			% ----------
			% A : Square, sparse transition matrix which does
            %	not include income or death transitions.
			%
			% g0 : Optional, the initial distribution for
			%	the iterative method.
			%
			% Returns
			% -------
			% g : The stationary distribution, of shape
			%	(nb_KFE, na_KFE, nz, ny).

			size_A = size(A);
			assert(size_A(1) == size_A(2), "KFESolver requires a square transition matrix")
			assert(numel(size_A) == 2, "KFESolver requires a square transition matrix")
			assert(issparse(A), "KFESolver requires a sparse transition matrix")

			if obj.options.iterate
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
			% Allows the user to manually set values in
			% the 'options' structure.
			if isprop(obj.options, keyword)
				obj.options.set(keyword, value);
			end
		end
	end

	methods (Access=private)

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
			while (iter <= obj.options.maxiters) && (dst > obj.options.tol)
				iter = iter + 1;

			    gg_tilde = obj.grdKFE.trapezoidal.diagm * g(:);
			    g1 = zeros(obj.p.nb_KFE*obj.p.na_KFE*obj.p.nz, obj.income.ny);
			    for iy = 1:obj.income.ny    
			    	gk_sum = sum(repmat(ytrans0p(iy,:), states_per_income, 1)...
			            .* reshape(gg_tilde, states_per_income, obj.income.ny),2);
			    	death_inflows = obj.compute_death_inflows(gg_tilde, iy);
		            g1(:,iy) = KFE_LHS{iy}*(gg_tilde(1+(iy-1)*states_per_income:iy*states_per_income)...
		                		+ obj.options.delta*gk_sum + obj.options.delta*death_inflows);
		        end

			    g1 = g1(:) ./ sum(g1(:));
			    g1 = obj.grdKFE.trapezoidal.diagm \ g1;

		        dst = max(abs(g1(:) - g(:)));

		        if obj.options.intermediate_check
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

				LHS{k} = (speye(states_per_income) - obj.options.delta * A(i1:i2, i1:i2)'...
			   		- obj.options.delta * (obj.income.ytrans(k,k) - obj.p.deathrate) * speye(states_per_income));
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
	    	if dst < obj.options.tol
			    fprintf('\tKFE converged after %i iterations\n', iter);
			elseif dst >= obj.options.tol
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