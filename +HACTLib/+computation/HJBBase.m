classdef (Abstract) HJBBase < handle

	properties (Abstract, Constant)
		required_parameters;
		required_income_vars;
	end

	properties (Constant)
		% Default options
		defaults = struct(...
				'implicit', false,...
				'delta', 1e5,...
				'HIS_maxiters', 0,...
				'HIS_tol', 1e-5,...
				'HIS_start', 2 ...
				);
 
	end

	properties (SetAccess=protected)
		%	An object with the following required attributes:
		%
		%		nb, na, nz
		%
		%		rho > 0
		%		- The time discount factor.
		%
		%		rhos
		%		- A vector indicating all possible values for rho.
		%		  If there is no rho heterogeneity, set rhos = [].
		%
		%		deathrate > 0
		p;

		%	An object with the following required attributes:
		%
		%		ny
		%	
		%		ytrans
		%		- The income transition matrix, should have row sums
		%		  of zero and shape (ny, ny).
		%
		%	and if using stochastic differential utility, the following
		%	is a required method:
		%
		%		income_transitions_SDU(p, V)
		%		- Must accept two arguments, p and V, and must return
		%		  an array of shape (nb*na*nz, ny, ny) to be used as
		%		  the adjusted income transition rates.
		%
		%	and--again for stochastic differential utility, the following
		%	method is required, or 'sparse_income_transitions' can be
		%	an attribute of the 'income' variable, in which case it must
		%	be a sparse array of size (nb*na*nz*ny, nb*na*nz*ny) indicating
		%	the adjusted income transition rates:
		%
		%		income_transition_matrix_SDU(p, ez_adj, 'HJB')
		%		- Must accept the three arguments above. The variable
		%		  'ez_adj' is the array generated by the
		%		  income_transitions_SDU() method.
		income;

		% Total number of states.
		n_states;

		% Number of states per income level.
		states_per_income;

		% Number of states along each dim, a vector.
		shape;

		% Sparse matrix of discount factors.
		rho_mat;
        
        % An HJBOptions object.
        options;

        % Current iteration, needed for the HIS.
        current_iteration = 0;
	end

	methods
		function obj = HJBBase(p, income, varargin)
			% See the properties block above or view this class'
			% documentation to see the requirements of the
			% input parameters.

			% ---------------------------------------------------------
			% Validate Input Arguments and Set Options
			% ---------------------------------------------------------
			obj.check_parameters(p);
			obj.p = p;

			obj.check_income(income);
			obj.income = income;

			obj.options = obj.parse_options(varargin{:});

			obj.n_states = p.nb * p.na * p.nz * income.ny;
			obj.states_per_income = p.nb * p.na * p.nz;
			obj.shape = [p.nb p.na p.nz income.ny];

			% ---------------------------------------------------------
			% Discount Factor Matrix
			% ---------------------------------------------------------
			obj.create_rho_matrix();
		end

		function V_update = solve(obj, A, u, V, varargin)
			% Updates the value function.

			obj.check_inputs(A, u, V);

			if obj.options.implicit
				V_update = obj.solve_implicit(A, u, V, varargin{:});
			else
				V_update = obj.solve_implicit_explicit(A, u, V, varargin{:});
			end
		end
	end

	methods (Access=protected)
		%%--------------------------------------------------------------
	    % Implicit Updating
	    % --------------------------------------------------------------
		function Vn1 = solve_implicit(obj, A, u, V, varargin)
			assert(obj.p.deathrate == 0,...
				"Fully implicit updating does not support death.")

			A_income = obj.income.full_income_transition_matrix(obj.p, V);

			if obj.returns_risk && obj.p.SDU
				risk_adj = reshape(vargin{1}, [], 1);
		        RHS = obj.options.delta * (u(:) + risk_adj) + Vn(:);
		    else
		    	RHS = obj.options.delta * u(:) + Vn(:);
		    end
	        
	        B = (obj.rho_mat - A - A_income) * obj.options.delta + speye(obj.n_states);
	        Vn1 = B \ RHS;
	        Vn1 = reshape(Vn1, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		%%--------------------------------------------------------------
	    % Implicit-Explicit Updating
	    % --------------------------------------------------------------
		function [Vn1, Bk_inv] = solve_implicit_explicit(obj, A, u, V, varargin)
			import HACTLib.computation.hjb_divisor
			obj.current_iteration = obj.current_iteration + 1;

			if obj.p.SDU && (obj.p.sigma_r > 0)
				sdu_adj = obj.income.income_transitions_SDU(obj.p, V);
				args1 = {sdu_adj};
				risk_adj_k = reshape(varargin{1}, [], obj.income.ny);
				args2 = {risk_adj_k};
			else
				args1 = {};
				args2 = {};
			end
			
			u_k = reshape(u, [], obj.income.ny);
			Vn_k = reshape(V, [], obj.income.ny);

			Vn1_k = zeros(obj.states_per_income, obj.income.ny);
			Bk_inv = cell(1, obj.income.ny);
			for k = 1:obj.income.ny

				[inctrans, inctrans_k] = obj.get_income_transitions(k, args1{:});

				Bk = hjb_divisor(obj.options.delta, obj.p.deathrate, k,...
					A, inctrans, obj.rho_mat);
            	Bk_inv{k} = inverse(Bk);

	        	Vn1_k(:,k) = obj.update_Vk_implicit_explicit(...
	        		Vn_k, u_k, k, Bk_inv{k}, inctrans_k, args2{:}, varargin{:});
	        end

	        Vn1 = reshape(Vn1_k, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		function Vn1_k = update_Vk_implicit_explicit(...
			obj, V_k, u_k, k, Bk_inv, inctrans_k, varargin)
        	
        	indx_k = ~ismember(1:obj.income.ny, k);

            offdiag_inc_term = sum(...
            	squeeze(inctrans_k(:,indx_k)) .* V_k(:,indx_k), 2);

            RHSk = obj.options.delta * (u_k(:,k) + offdiag_inc_term)...
            		+ V_k(:,k);
            
            if obj.p.SDU && (obj.p.sigma_r > 0)
            	risk_adj = varargin{1};
           		RHSk = RHSk + obj.options.delta * risk_adj(:,k);
           	end
            
            Vn1_k = Bk_inv * RHSk;
        end

		function options = parse_options(obj, varargin)
			import HACTLib.computation.HJBBase
			import HACTLib.aux.parse_keyvalue_pairs
			import HACTLib.Checks

			defaults = HJBBase.defaults;
			options = parse_keyvalue_pairs(defaults, varargin{:});

			mustBePositive(options.delta);
			mustBePositive(options.HIS_tol);

			super_class = split(class(obj), '.');
		    super_class = super_class(end);

			Checks.is_integer(super_class, options.HIS_maxiters);
			Checks.is_logical(super_class, options.implicit);
		end

		function obj = create_rho_matrix(obj)
			if obj.options.implicit
				% discount factor values
		        if numel(obj.p.rhos) > 1
		            rhocol = repmat(kron(obj.p.rhos(:), ones(obj.p.nb*obj.p.na, 1)), obj.income.ny, 1);
		            obj.rho_mat = spdiags(rhocol, 0, obj.n_states, obj.n_states);
		        else
		            obj.rho_mat = obj.p.rho * speye(obj.n_states);
		        end
		    else
		    	if numel(obj.p.rhos) > 1
			        rhocol = kron(obj.p.rhos(:), ones(obj.p.nb*obj.p.na,1));
			        obj.rho_mat = spdiags(rhocol, obj.states_per_income, obj.states_per_income);
			    else
			        obj.rho_mat = obj.p.rho * speye(obj.states_per_income);
			    end
	    	end
		end

		function check_inputs(obj, A, u, V)
			import HACTLib.Checks;
            
            super_class = split(class(obj), '.');
            super_class = super_class(end);

			Checks.is_square_sparse_matrix(super_class, A, obj.n_states);
			Checks.has_shape(super_class, u, obj.shape);
			Checks.has_shape(super_class, V, obj.shape);
		end
	end

	methods (Abstract, Access=protected)
		check_if_SDU(obj);
		[inctrans, inctrans_k] = get_income_transitions(obj, varargin);
		check_parameters(obj, p);
		check_income(obj, income);
	end
end