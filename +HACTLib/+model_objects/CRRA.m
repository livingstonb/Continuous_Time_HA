classdef CRRA
	% Provides functions for CRRA utility, with or
	% without IES heterogeneity.

	methods (Static)
		function u = utility(c, invies)
			import HACTLib.aux.replace_where

			u_log = log(c);
			u_other = c .^ (1 - invies) ./ (1 - invies);
			u = replace_where(invies==1, u_log, u_other);
		end

		function muc = marginal_utility(c, invies)
			muc = c .^ (-invies); 
		end

		function c = u1inv(v, invies, zdim)
			c = v .^ (-1 ./ invies);
		end
	end
end