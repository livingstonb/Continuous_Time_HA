classdef HJBSolverSDU < HACT_Tools.algorithms.HJBBase

	methods
		function obj = HJBSolverSDU(p, income, options)
			obj = obj@HACT_Tools.algorithms.HJBBase(p, income, options);
		end
	end

	properties (Constant)
		% Update this array when the required parameters
		% change.
		required_parameters = {'nb', 'na', 'nz',...
			'rho', 'rhos', 'deathrate'};

		% Update this array when the required income
		% variables change.
		required_income_vars = {'ny', 'ytrans'};
	end

	methods (Access=protected)
		function check_if_SDU(obj)
			msg = [	"Model does not use stochastic differential utility, ",...
					"this subclass should not be used."]
			assert(obj.p.SDU, msg);
		end

		%%--------------------------------------------------------------
	    % Implicit Updating
	    % --------------------------------------------------------------
		function Vn1 = solve_implicit(obj, A, u, V, varargin)
			assert(obj.p.deathrate == 0,...
				"Fully implicit updating does not support death.")

			A_income = obj.income.income_transition_matrix_SDU(obj.p, V);

			if obj.returns_risk
				risk_adj = reshape(vargin{1}, [], 1);
		        RHS = obj.options.delta * (u(:) + risk_adj) + Vn(:);
		    else
		    	RHS = obj.options.delta * u(:) + Vn(:);
		    end
	        
	        B = (obj.rho_mat - A - A_income) * obj.options.delta + speye(obj.n_states);
	        Vn1 = B \ RHS;
	        Vn1 = reshape(Vn1, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		function Vn1 = solve_implicit_explicit(obj, A, u, V, varargin)
			obj.current_iteration = obj.current_iteration + 1;

			u_k = reshape(u, [], obj.income.ny);
			Vn_k = reshape(V, [], obj.income.ny);
			risk_adj_k = reshape(varargin{1}, [], obj.income.ny);

			ez_adj = obj.income.income_transitions_SDU(obj.p, V);

			Vn1_k = zeros(obj.states_per_income, obj.income.ny);			
			for k = 1:obj.income.ny
				inctrans = spdiags(ez_adj(:, k, k), 0,...
					obj.states_per_income, obj.states_per_income);
				Bk = obj.construct_Bk(k, A, inctrans, varargin{:});
            	Bk_inv = inverse(Bk);
	        	Vn1_k(:,k) = obj.update_Vk_implicit_explicit(...
	        		Vn_k, u_k, k, Bk_inv, ez_adj, risk_adj_k(:,k), varargin{:});
	        end
	        Vn1 = reshape(Vn1_k, obj.shape);
		end

		function Vn1_k = update_Vk_implicit_explicit(...
			obj, V_k, u_k, k, Bk_inv, ez_adj, varargin)
        	
        	indx_k = ~ismember(1:obj.income.ny, k);

            offdiag_inc_term = sum(...
            	squeeze(ez_adj(:, k, indx_k)) .* V_k(:, indx_k), 2);

            RHSk = obj.options.delta * u_k(:,k)...
            	+ V_k(:,k) + obj.options.delta * offdiag_inc_term;
           	
           	if numel(varargin) == 1
           		RHSk = RHSk + obj.options.delta * varargin{1};
           	end
            
            Vn1_k = Bk_inv * RHSk;
        end
	end
end