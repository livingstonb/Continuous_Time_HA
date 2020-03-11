classdef AdjustmentCost
	properties
		a_lb;
		kappa0;
		kappa1;
		kappa2;
	end

	methods
		function set_form1(obj, a_lb, chi0, chi1, chi2)
			obj.a_lb = a_lb;
			obj.kappa0 = chi0;
			obj.kappa1 = chi1 ^ (-chi2);
			obj.kappa2 = chi2;
		end

		function set_form2(obj, a_lb, kappa0, kappa1, kappa2)
			obj.a_lb = a_lb;
			obj.kappa0 = kappa0;
			obj.kappa1 = kappa1;
			obj.kappa2 = kappa2;
		end

		function cost_out = compute_cost(d, a_grid)
			a_scaled = max(a_grid, obj.a_lb);

			cost_linear = obj.kappa0 * abs(d);
			cost_concave = obj.kappa1 / (1 + obj.kappa2) ...
				* abs(d) .^ (1 + obj.kappa2) ...
				./ (a_scaled .^ obj.kappa2);

			cost_out = cost_linear + cost_concave;
		end

		% function cost_deriv = compute_deriv(d, a_grid)
		% 	deriv_dneg = - obj.kappa0 - 
		% end

		function d_opt = opt_deposits(Vb, Va, a)
			dpos_term = Va ./ Vb - 1 - obj.kappa;
			dneg_term = 1 - obj.kappa0 - Va ./ Vb;

			a_scaled = max(a, obj.a_lb);
			d_opt = zeros(size(a));
			d_opt(dpos_term >0) = a_scaled .* dpos_term .^ (1/obj.kappa2);
			d_opt(dneg_term > 0) = - a_scaled .* dneg_term .^ (1/obj.kappa2);
		end
	end

	methods (Access=protected)

		function cost_out = form1_cost_handle(d, a_grid, a_lb,...
			chi0, chi1, chi2)
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

			linear_term = chi0 * abs(d);

			a_scaled = max(a, a_lb);
			coeff = (chi1 ^ (-chi2)) / (1 + chi2);
			concave_term = coeff 
			% concave_term = (chi1 / (1 + chi2)) * abs(d) .^ (1 + chi2) ...
			% 	.* a_scaled .^ chi2;

			% cost_out = linear_term + concave_term;


			d_scaled = d./max(a_grid,p.a_lb);
    		adj_cost = max(a_grid,p.a_lb) .* (p.chi0 * abs(d_scaled) + 1/(1+p.chi2) * (abs(d_scaled).^(1+p.chi2) * p.chi1^(-p.chi2)));
		end

		function cost_out = form2_cost_handle(d, a, a_lb,...
			kappa0, kappa1, kappa2)
			% Parameters
			% ----------
			% d : deposit rate
			%
			% a_grid : illiquid asset levels
			%
			% Returns
			% -------
			% adj_cost : the adjustment cost, same shape as d

			linear_term = kappa0 * abs(d);

			a_scaled = max(a, a_lb);
			concave_term = (kappa1 / (1 + kappa2)) * abs(d) .^ (1 + kappa2) ...
				./ (a_scaled .^ kappa2);

			cost_out = linear_term + concave_term;	
		end
	end

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

			linear_term = sign(d) * p.chi0;
			power_term = sign(d) .* (abs(d_scaled) / p.chi1) .^ p.chi2;
			chi_prime = linear_term + power_term;
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