classdef AdjustmentCost
	% static methods for the adjustment cost function and related functions.

	methods(Static)
		function adj_cost = cost(d, a_grid, p)
			% adjustment cost function chi(d)

			% Parameters
			% ----------
			% d : deposit rate
			%
			% a_grid : illiquid asset levels
			%
			% p : a Params object
			%
			% Returns
			% -------
			% adj_cost : the adjustment cost, same shape as d

			d_scaled = d./max(a_grid,p.a_lb);
    		adj_cost = max(a_grid,p.a_lb) .* (p.chi0 * abs(d_scaled) + 1/(1+p.chi2) * (abs(d_scaled).^(1+p.chi2) * p.chi1^(-p.chi2)));
		end

		function chi_prime = derivative(d, a_grid, p)
			% derivative of the adjustment cost function, chi'(d)

			% Parameters
			% ----------
			% d : deposit rate
			%
			% a_grid : illiquid asset levels
			%
			% p : a Params object
			%
			% Returns
			% -------
			% chi_prime : the derivative, same shape as d

			d_scaled = d ./ max(a_grid, p.a_lb);
			chi_prime = p.chi0 + sign(d) .* p.chi1 ^(-p.chi2) .* abs(d_scaled) .^ p.chi2;
		end

		function d = derivative_inverse(y, a_grid, p)
		    % inverse of the derivative of the adjustment cost function,
		    % i.e. (chi'(d))^{-1}

		    % Parameters
			% ----------
			% y : input values at which to evaluate the functikon
			%
			% a_grid : illiquid asset levels
			%
			% p : a Params object
			%
			% Returns
			% -------
			% d : the inverse of the derivative of the adjustment cost function,
			%	same shape as y

			d = sign(y) .* max(a_grid,p.a_lb) .* p.chi1 .* (abs(y) - p.chi0).^(1/p.chi2);
		end
	end
end