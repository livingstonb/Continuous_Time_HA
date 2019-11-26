classdef HJBSolverSDU < HACT_Tools.algorithms.HJBSolver
	methods (Access=protected)
		function check_if_SDU(obj)
			assert(obj.p.SDU,...
				[	"Model does not use stochastic differential utility, ",...
					"this subclass should not be used."])
		end

		function RHS = get_RHS_implicit(obj, u, V, varargin)
			% Overloads the same method in HJBSolver.
			% The risk adjustment term should be passed as
			% a third argument.

			if nargin == 3
				RHS = get_RHS_implicit@HACT_Tools.algorithms.HJBSolver(obj, u, V);
			elseif nargin == 4
				HACT_Tools.Checks.have_same_shape(V, varargin{1});
				RHS = obj.options.delta * (u(:) + reshape(varargin{1}, [], 1)) + Vn(:);
			else
				error("Expected <= 3 arguments");
			end
		end

		function [Vn1_k, Bk_inv] = loop_over_income_implicit_explicit(obj, V, A, u, varargin)
			narginchk(4, 5);
        	u_k = reshape(u, [], obj.income.ny);
        	Vn_k = reshape(V, [], obj.income.ny);

        	ez_adj = obj.income.income_transitions_SDU(obj.p, V);

        	if nargin > 4
        		risk_adj_k = reshape(varargin{1}, [], obj.income.ny);
        	end

        	Vn1_k = zeros(obj.states_per_income, obj.income.ny);
        	Bk_inv = cell(1, obj.income.ny);
        	for k = 1:obj.income.ny
        		Bk = obj.construct_Bk(k, A, ez_adj(:, k, k));
            	Bk_inv{k} = inverse(Bk);

            	indx_k = ~ismember(1:obj.income.ny, k);

	            offdiag_inc_term = sum(squeeze(ez_adj(:, k, indx_k))...
	           		.* Vn_k(:, indx_k), 2);

	            RHSk = obj.options.delta * u_k(:,k)...
		            	+ Vn_k(:,k) + obj.options.delta * offdiag_inc_term;
		            	
            	if nargin > 4
            		RHSk = RHSk + obj.options.delta * risk_adj_k(:, k);
            	end
	            
	            Vn1_k(:,k) = Bk_inv{k} * RHSk;
        	end
	    end

	    function inc_term = get_Bk_income_term(obj, k, varargin)
	    	% Adds risk-adjusted diagonal income transitions to
	    	% the left-hand-side of the HJB.

	    	inc_term = spdiags(varargin{1}, 0,...
	    		obj.states_per_income,...
	    		obj.states_per_income);
	    end
	end
end