classdef Preferences < handle

	properties
		% Utility function
		u;

		% Marginal utility
		u1;

		% Inverse of marginal utility function
		u1inv;

		% Labor disutility
		hrs_u;

		% Derivative of labor disutility
		hrs_u1;

		% Inverse of marginal labor disutility
		hrs_u1inv;
	end

	methods
		function set_crra(obj, riskaver)
			import HACTLib.model_objects.CRRA.*
			obj.u = @(c) HACT
		end

		function set_crra_het(obj, riskavers)
		end
	end
end