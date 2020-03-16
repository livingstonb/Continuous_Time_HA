classdef KernelSmoother < handle
	properties
		x;
		x_hat;
		y;
		h;
	end

	methods
		function obj = KernelSmoother(x, y, varargin)
			if nargin > 0
				obj.x = x;
				obj.y = y;
				obj.set_bandwidth(varargin{:});
			end
		end

		function set_from_pmf(obj, pmf, values, varargin)
			arr_sorted = sortrows([values(:), pmf(:)]);

			[v_u, iu] = unique(arr_sorted(:,1), 'last');
			cdf_u = cumsum(arr_sorted(:,2));
			cdf_u = cdf_u(iu);

			obj.x = v_u;
			obj.y = cdf_u;

			obj.set_bandwidth(varargin{:});
		end

		function x_out = keval(obj, x_query)
			x_query = reshape(x_query, 1, []);
			kd = exp(-(x_query - obj.x) .^ 2 ./ (2 * obj.h .^ 2));
			
			x_out = sum(obj.y .* kd, 1) ./ sum(kd, 1);
			x_out = x_out(:);
        end
        
   %      function x_out = keval_epanechnikov(obj, x_query)
			% x_query = reshape(x_query, 1, []);

			% v = abs(x_query - obj.x) ./ obj.h;
			% kd = 1 - v .^ 2;
			% kd(v >= 1) = 0;
			
			% x_out = sum(obj.y .* kd, 1) ./ sum(kd, 1);
			% x_out = x_out(:);
   %      end

        function y_hat = keval_inv(obj, y_query)
        	y_grid = obj.keval(obj.x);

        	ginterp = griddedInterpolant(y_grid(obj.x>0), obj.x(obj.x>0),...
        		'spline', 'nearest');
        	y_hat = ginterp(y_query);
        end

        function make_plots(obj)
        	plot(obj.x, obj.y)
            hold on

            plot(obj.x, obj.keval(obj.x));
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