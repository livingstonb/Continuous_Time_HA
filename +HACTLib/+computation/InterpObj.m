
classdef InterpObj < handle
	properties
		x;
		y;

		kernel_smoothing = false;
	end

	properties (Access=protected)
		kernel_smoother;
		pct_interp;
		cdf_interp;
	end

	methods
		function obj = InterpObj()
		end

		function set_dist(obj, values, pmf, dims_to_keep)
			import HACTLib.aux.multi_sum

			if nargin < 4
				already_sorted = true;
			else
				already_sorted = isequal(dims_to_keep, []);
			end

			if already_sorted
				[v_u, iu] = unique(values(:), 'last');
				cdf_u = cumsum(pmf(:));
			else
				sum_dims = 1:4;
				sum_dims = sum_dims(...
					~ismember(sum_dims, dims_to_keep));

				pmf_x = multi_sum(pmf, sum_dims);

				sorted_mat = sortrows([values(:), pmf_x(:)]);
				[v_u, iu] = unique(sorted_mat(:,1), 'last');
				cdf_u = cumsum(sorted_mat(:,2));
			end

			if numel(v_u) > 5000
				[obj.x, obj.y] = thin_cdf(v_u, cdf_u(iu), 5000);
			else
				obj.x = v_u;
				obj.y = cdf_u(iu);
			end
		end

		function configure(obj, kernel_options)
			if numel(kernel_options) == 0
				obj.configure_interpolants();
			else
				obj.kernel_smoothing = true;
				obj.configure_kernel_smoother(kernel_options);

				obj.x = [];
				obj.y = [];
			end
		end

		function configure_interpolants(obj)
			% Interpolant for cdf
			obj.cdf_interp = griddedInterpolant(...
				obj.x, obj.y, 'pchip', 'nearest');

			% Interpolant for percentiles
			[y_u, iu] = unique(obj.y);
			x_u = obj.x(iu);
			obj.pct_interp = griddedInterpolant(...
				y_u, x_u, 'pchip', 'nearest');
		end

		function configure_kernel_smoother(obj, kernel_options)
			import HACTLib.computation.KernelSmoother

			obj.kernel_smoother = KernelSmoother(kernel_options);
			obj.kernel_smoother.set(obj.x, obj.y);
		end


		function z_out = icdf(obj, p)
			if obj.kernel_smoothing
				z_out = obj.kernel_smoother.keval_inv(p);
			else
				z_out = obj.pct_interp(p);
			end
		end

		function cdf_out = cdf(obj, val)
			if obj.kernel_smoothing
				cdf_out = obj.kernel_smoother.keval(val);
			else
				cdf_out = obj.cdf_interp(val);
			end
		end

		function plot_cdf(obj, varargin)
			if obj.kernel_smoothing
				obj.kernel_smoother.make_plots(varargin{:});
			else
				plot(obj.x, obj.y);
			end
			ylim([0, 0.8]);
			xlim('auto')
		end
	end
end

function [x, y] = thin_cdf(x, y, nmax)
	n = numel(x);

	pct_to_drop = 100 * (n - nmax) / n;

	% xn = (x - min(x)) / (max(x) - min(x));
	% xn = log(0.01 + xn);
	% yn = 100 * (y - min(y)) / (max(y) - min(y));

	dydx = diff(y) ./ diff(x);
	dydx_pctile = prctile(dydx, pct_to_drop);
	keep = (dydx > dydx_pctile);

	% dnorm = sqrt(sum(diff([xn, yn]) .^ 2, 2));
	% dnorm = [100; dnorm];

	% dx = prctile(dnorm, pct_to_drop);

	% keep = (dnorm >= dx);
    x = x(keep);
    y = y(keep);
end