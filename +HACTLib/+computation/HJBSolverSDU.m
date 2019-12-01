classdef HJBSolverSDU < HACTLib.computation.HJBBase

	methods
		function obj = HJBSolverSDU(p, income, options)
			obj = obj@HACTLib.computation.HJBBase(p, income, options);
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
					"this subclass should not be used."];
			assert(obj.p.SDU, msg);
		end

		function [inctrans, inctrans_k] = get_income_transitions(obj, k, sdu_adj, varargin)
			inctrans = HACTLib.aux.sparse_diags(sdu_adj(:,k,k), 0);
			inctrans_k = squeeze(sdu_adj(:,k,:));
		end

        function check_parameters(obj, p)
			HACTLib.Checks.has_attributes(...
				'HJBSolverSDU', p, obj.required_parameters);
		end

		function check_income(obj, income)
			HACTLib.Checks.has_attributes('HJBSolverSDU',...
				income, obj.required_income_vars);
			HACTLib.Checks.is_square_matrix('HJBSolverSDU',...
				income.ytrans);
		end
	end
end