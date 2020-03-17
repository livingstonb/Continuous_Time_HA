classdef KernelSmoother < handle
	properties
		x;
		y_hat;
		y;
		h;
		ktype = 'gaussian';
		log_transform = false;
		transform_x = @(x) x;
	end

	methods
		function obj = KernelSmoother(ktype, log_transform)
			if nargin == 1
				obj.ktype = ktype;
			elseif nargin == 2
				obj.ktype = ktype;
				obj.log_transform = log_transform;

				if log_transform
					obj.transform_x = @(x) exp(x) - 0.1;
				end
			end
		end

		function set_for_cdf_smoothing(obj, pmf, values, h)
			arr_sorted = sortrows([values(:), pmf(:)]);

			[v_u, iu] = unique(arr_sorted(:,1), 'last');
			cdf_u = cumsum(arr_sorted(:,2));
			cdf_u = cdf_u(iu);

			if obj.log_transform
				v_u = log(0.1 + v_u);
			end

			obj.x = v_u;
			obj.y = cdf_u;

			obj.h = h;
			obj.y_hat = impose_monotonic(obj.x, obj.keval(obj.x));
		end

		% function set_for_pmf_smoothing(obj, pmf, values, h)
		% 	obj.x = v_u;
		% 	obj.y = arr_sorted(iu,2);

		% 	obj.h = h;
		% 	obj.y_hat = impose_monotonic(obj.x, obj.keval(obj.x));
		% end

		function set(obj, x, y, h)
			if obj.log_transform
				obj.x = log(0.1 + x);
			else
				obj.x = x;
			end
			obj.y = y;
			obj.h = h;
		end

		function x_out = keval(obj, x_query)
			x_query = reshape(x_query, 1, []);
			v = abs(x_query - obj.x) ./ obj.h;

			switch obj.ktype
				case 'gaussian'
					kd = exp(-v .^ 2);
				case 'epanechnikov'
					kd = 1 - v .^ 2;
					kd(v > 1) = 0;
				case 'logistic'
					kd = 1 ./ (exp(v) + 2 + exp(-v));
				case 'triweight'
					kd = (1 - v .^ 2) .^ 3;
					kd(v > 1) = 0;
				end
			
			x_out = sum(obj.y .* kd, 1) ./ sum(kd, 1);
			x_out = x_out(:);
        end
        
        function x_hat = keval_inv(obj, y_query)
        	ginterp = griddedInterpolant(obj.y_hat, exp(obj.x),...
        		'spline', 'nearest');
        	x_hat = ginterp(y_query);
        end

        function make_plots(obj)
        	if obj.log_transform
        		x_grid = obj.transform_x(obj.x);
        	else
        		x_grid = obj.x;
        	end

	        plot(x_grid, obj.y)
	        hold on
            plot(x_grid, obj.keval(obj.x));
            legend('Raw data', 'Kernel-smoothed (Gaussian)')
        end
        
        function plot_cdfs(obj)
            plot(obj.x, cumsum(obj.y))
            hold on

            plot(obj.x, cumsum(obj.keval(obj.x)));
            legend('Raw data', 'Kernel-smoothed (Gaussian)')
        end
	end

	methods (Access=protected)
		function set_bandwidth(obj, varargin)
			if nargin == 1
				obj.h = sqrt(obj.x) / max(obj.x);
				% obj.h = max(obj.h, min(obj.h(obj.h>0)));
			else
				obj.h = varargin{1};
			end
		end
	end
end

function arr_out = impose_monotonic(agrid, arr_in)
	increasing = false(size(arr_in));

	last_val = arr_in(1);
	increasing(1) = true;
	for ia = 2:numel(arr_in)
		if arr_in(ia) > last_val
			increasing(ia) = true;
			last_val = arr_in(ia);
		end
	end

	cinterp = griddedInterpolant(agrid(increasing),...
		arr_in(increasing), 'spline', 'nearest');

	arr_out = zeros(size(arr_in));
	arr_out(~increasing) = cinterp(agrid(~increasing));
	arr_out(increasing) = arr_in(increasing);
end