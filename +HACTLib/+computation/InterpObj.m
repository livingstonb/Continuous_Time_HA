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

			sum_dims = 1:4;
			sum_dims = sum_dims(...
				~ismember(sum_dims, dims_to_keep));

			pmf_x = multi_sum(pmf, sum_dims);

			sorted_mat = sortrows([values(:), pmf_x(:)]);
			[v_u, iu] = unique(sorted_mat(:,1), 'last');
			cdf_u = cumsum(sorted_mat(:,2));

			obj.x = v_u;
			obj.y = cdf_u(iu);
		end

		function configure(obj, kernel_options)
			if nargin == 1
				obj.configure_interpolants();
			else
				obj.kernel_smoothing = true;
				obj.configure_kernel_smoother(kernel_options);
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

			log_transform = true;
			obj.kernel_smoother = KernelSmoother(...
				kernel_options.ktype, kernel_options.log_transform);
			obj.kernel_smoother.set(obj.x, obj.y, kernel_options.h);
		end


		function val_pct = percentile(obj, pct)
			if obj.kernel_smoothing
				val_pct = obj.kernel_smoother.keval_inv(pct);
			else
				val_pct = obj.pct_interp(pct);
			end
		end

		function cdf_out = cdf_at_value(obj, val)
			if obj.kernel_smoothing
				cdf_out = obj.kernel_smoother.keval(val);
			else
				cdf_out = obj.cdf_interp(val);
			end
		end

		function plot_cdf(obj)
			plot(obj.x, obj.y);
			ylim([0, 0.8]);
			xlim('auto')
		end
	end
end