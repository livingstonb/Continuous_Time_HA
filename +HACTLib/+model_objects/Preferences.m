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
		function set_crra(obj, invies)
			import HACTLib.model_objects.CRRA
			narginchk(2, 3);

			obj.u = @(c) CRRA.utility(c, invies);
			obj.u1 = @(c) CRRA.marginal_utility(c, invies);
			obj.u1inv = @(u) CRRA.u1inv(u, invies);
		end

		function set_sdu(obj, invies)
			obj.u = @(c, rho) rho .* CRRA.utility(c, invies);
			obj.u1 = @(c, rho) rho .* CRRA.marginal_utility(c, invies);
			obj.u1inv = @(u, rho) CRRA.u1inv(u ./ rho, invies);
		end
	end
end