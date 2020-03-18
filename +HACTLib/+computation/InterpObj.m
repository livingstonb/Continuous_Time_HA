
classdef InterpObj < handle
	properties (Access=protected)
		x;
		y;
		kernel_smoother;
		linear_interpolant;
		insufficient_grid_pts = false;
		kernel_smoothing = false;
	end

	methods
		function obj = InterpObj()
		end

		function set_dist(obj, values, pmf, dims_to_keep)
			import HACTLib.aux.multi_sum

			if (nargin == 3) || isequal(dims_to_keep, []);
				% Already sorted
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
			elseif numel(v_u) <= 5
				obj.insufficient_grid_pts = true;
			else
				obj.x = v_u;
				obj.y = cdf_u(iu);
			end
		end

		function configure(obj, kernel_options)
			if obj.insufficient_grid_pts
				return
			end

			obj.configure_interpolants();
			
			if numel(kernel_options) == 1
				obj.kernel_smoothing = true;
				obj.configure_kernel_smoother(kernel_options);
			end
		end

		function configure_interpolants(obj)
			obj.linear_interpolant = struct();

			% Interpolant for cdf
			obj.linear_interpolant.cdf = griddedInterpolant(...
				obj.x, obj.y, 'pchip', 'nearest');

			% Interpolant for inv cdf
			[y_u, iu] = unique(obj.y);
			x_u = obj.x(iu);
			obj.linear_interpolant.icdf = griddedInterpolant(...
				y_u, x_u, 'pchip', 'nearest');
		end

		function configure_kernel_smoother(obj, kernel_options)
			import HACTLib.computation.KernelSmoother

			obj.kernel_smoother = KernelSmoother(kernel_options);
			obj.kernel_smoother.set(obj.x, obj.y);
		end


		function z_out = icdf(obj, p, interp_type)
			if obj.insufficient_grid_pts
				z_out = NaN;
			elseif (nargin == 2) && (obj.kernel_smoothing)
				z_out = obj.kernel_smoother.keval_inv(p);
			elseif (nargin == 2)
				z_out = obj.linear_interpolant.icdf(p);
			elseif (nargin == 3) && strcmp(interp_type, 'kernel')
				z_out = obj.kernel_smoother.keval_inv(p);
			elseif (nargin == 3) && strcmp(interp_type, 'linear')
				z_out = obj.linear_interpolant.icdf(p);
			end
		end

		function g_out = cdf(obj, val, interp_type)
			if obj.insufficient_grid_pts
				g_out = NaN;
			elseif (nargin == 2) && (obj.kernel_smoothing)
				g_out = obj.kernel_smoother.keval(val);
			elseif (nargin == 2)
				g_out = obj.linear_interpolant.cdf(val);
			elseif (nargin == 3) && strcmp(interp_type, 'kernel')
				g_out = obj.kernel_smoother.keval(val);
			elseif (nargin == 3) && strcmp(interp_type, 'linear')
				g_out = obj.linear_interpolant.cdf(val);
			end
		end

		function g_out = cdf_spline(obj, val, show_plot)
			if obj.insufficient_grid_pts
				g_out = NaN;
			end

			if nargin == 2
				show_plot = false;
			end

			lsp = HACTLib.computation.LocalSpline();
			g_out = lsp.fit_spline(obj.x, obj.y, val);

			if show_plot
				lsp.show_fit(obj.x, obj.y);
			end
		end

		% function g_out = cdf_poly(obj, val, interval)
		% 	keep = (obj.x >= interval(1)) & ...
		% 		(obj.x <= interval(2));
		% 	p = polyfit(obj.x(keep), obj.y(keep), 3);
		% 	g_out = polyval(p, val);
		% end

		function plot_cdf(obj, varargin)
			if obj.insufficient_grid_pts
				error('Insufficient number of grid points')
			elseif obj.kernel_smoothing
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

	dydx = diff(y) ./ diff(x);
	dydx_pctile = prctile(dydx, pct_to_drop);
	keep = (dydx > dydx_pctile);
	keep = [true; keep];

    x = x(keep);
    y = y(keep);
end