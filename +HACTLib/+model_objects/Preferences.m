classdef Preferences < handle

	properties
		% Utility, u(c)
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

		% Utility, u(c, h)
		u_total;
	end

	methods
		function set_crra(obj, invies)
			import HACTLib.model_objects.CRRA

			obj.u = @(c) CRRA.utility(c, invies);
			obj.u1 = @(c) CRRA.marginal_utility(c, invies);
			obj.u1inv = @(u) CRRA.u1inv(u, invies);
		end

		function set_frisch(obj, coeff, frisch)
			import HACTLib.model_objects.Frisch

			obj.hrs_u = @(h) Frisch.labor_disutility(h, coeff, frisch);
			obj.hrs_u1 = @(h) Frisch.marginal_disutility(h, coeff, frisch);
			obj.hrs_u1inv = @(h) Frisch.inv_marginal_disutility(h, coeff, frisch);
		end
	end
end