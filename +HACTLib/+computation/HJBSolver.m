classdef HJBSolver < HACTLib.computation.HJBBase
	% A solver for the Hamilton-Jacobi-Bellman equation.
	% Implicit and implicit-explicit solution methods are
	% provided.
	%
	% The recommended use of this class is to first create
	% a Params object and an Income object, and pass these
	% to the HJBSolver constructor. Alternatively, one can
	% use other objects that satisfy the requirements laid
	% out in the properties block of HJBBase.

	properties (Constant)
		% Update this array when the required parameters
		% change.
		required_parameters = {'nb', 'na', 'nz',...
			'rho', 'rhos', 'deathrate'};

		% Update this array when the required income
		% variables change.
		required_income_vars = {'ny', 'ytrans'};
	end

	methods
		function obj = HJBSolver(p, income, options)
			obj = obj@HACTLib.computation.HJBBase(p, income, options);
		end
	end

	methods (Access=protected)
		function check_if_SDU(obj)
			msg = [	"Model uses stochastic differential utility, "...
					"the appropriate subclass should be used"];
			assert(~obj.p.SDU, msg);
		end

		%%--------------------------------------------------------------
	    % Implicit Updating
	    % --------------------------------------------------------------
		function Vn1 = solve_implicit(obj, A, u, V, varargin)
			assert(obj.p.deathrate == 0,...
				"Fully implicit updating does not support death.")

	        % Add income transitions
	        B = A + kron(obj.income.ytrans, speye(obj.states_per_income));

	        RHS = obj.options.delta * u(:) + Vn(:);
	        
	        B = (obj.rho_mat - B) * obj.options.delta + speye(obj.n_states);
	        Vn1 = B \ RHS;
	        Vn1 = reshape(Vn1, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		%%--------------------------------------------------------------
	    % Implicit-Explicit Updating
	    % --------------------------------------------------------------
		function Vn1 = solve_implicit_explicit(obj, A, u, V, varargin)
			import HACTLib.computation.hjb_divisor
			obj.current_iteration = obj.current_iteration + 1;

			u_k = reshape(u, [], obj.income.ny);
			Vn_k = reshape(V, [], obj.income.ny);

			Vn1_k = zeros(obj.states_per_income, obj.income.ny);
			Bk_inv = cell(1, obj.income.ny);
			for k = 1:obj.income.ny

				inctrans = obj.income.ytrans(k, k) * speye(obj.states_per_income);
				Bk = hjb_divisor(obj.options.delta, obj.p.deathrate, k,...
					A, inctrans, obj.rho_mat);
            	Bk_inv{k} = inverse(Bk);
	        	Vn1_k(:,k) = obj.update_Vk_implicit_explicit(...
	        		Vn_k, u_k, k, Bk_inv{k}, varargin{:});
	        end

	        % Howard improvement step
	        if obj.current_iteration >= obj.options.HIS_start
	            Vn1_k = obj.howard_improvement_step(Vn1_k, u_k, Bk_inv);
	        end
	        Vn1 = reshape(Vn1_k, obj.p.nb, obj.p.na, obj.p.nz, obj.income.ny);
		end

		function Vn1_k = update_Vk_implicit_explicit(...
			obj, V_k, u_k, k, Bk_inv, varargin)
        	
        	indx_k = ~ismember(1:obj.income.ny, k);

            offdiag_inc_term = sum(...
            	repmat(obj.income.ytrans(k,indx_k), obj.states_per_income, 1)...
                .* V_k(:,indx_k), 2);

            RHSk = obj.options.delta * u_k(:,k)...
            	+ V_k(:,k) + obj.options.delta * offdiag_inc_term;
            
            Vn1_k = Bk_inv * RHSk;
	    end

		function Vn2_k = howard_improvement_step(obj, Vn1_k, u_k, Bik_all)
		    % Technique to speed up convergence.

		    for jj = 1:obj.options.HIS_maxiters
		        Vn2_k = zeros(obj.states_per_income, obj.income.ny);
		        for kk = 1:obj.income.ny
		            indx_k = ~ismember(1:obj.income.ny, kk);
		            
		            Vkp_stacked = sum(...
		                	repmat(obj.income.ytrans(kk, indx_k), obj.states_per_income,1)...
		                	.* Vn1_k(:,indx_k), 2);
		            qk = obj.options.delta * (u_k(:,kk) + Vkp_stacked) + Vn1_k(:,kk);
		            Vn2_k(:,kk) = Bik_all{kk} * qk;
		        end

		        dst = max(abs(Vn2_k(:) - Vn1_k(:)));
		        Vn1_k = Vn2_k;
		        if dst < obj.options.HIS_tol
		            break
		        end
    		end
		end

		% function check_inputs(obj, A, u, V)
		% 	check_inputs@HACTLib.computation.HJBBase(class(obj), A, u, V);
		% end

		function check_parameters(obj, p)
			HACTLib.Checks.has_attributes(...
				'HJBSolver', p, obj.required_parameters);
		end

		function check_income(obj, income)
			HACTLib.Checks.has_attributes('HJBSolver',...
				income, obj.required_income_vars);
			HACTLib.Checks.is_square_matrix('HJBSolver',...
				income.ytrans);
		end
	end
end