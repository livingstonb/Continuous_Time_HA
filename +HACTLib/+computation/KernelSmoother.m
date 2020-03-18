classdef KernelSmoother < handle
	properties (Constant)
		defaults = struct(...
			'ktype', 'gaussian',...
			'force_fit_cdf_low', [],...
			'rescale_and_log', true,...
			'h', 0.2...
			)
	end

	properties
		x;
		inv_interp;
		y;
		h;

		ktype;
        rescale_and_log;
		force_fit_cdf_low;
		x_transform;

		p_lower;
	end

	methods
		function obj = KernelSmoother(varargin)
			import HACTLib.aux.parse_keyvalue_pairs
			options = parse_keyvalue_pairs(obj.defaults, varargin{:});

			obj.ktype = options.ktype;
			obj.force_fit_cdf_low = options.force_fit_cdf_low;
            obj.rescale_and_log = options.rescale_and_log;
			obj.h = options.h;
		end

		function set(obj, x, y)
			obj.x = x;
			obj.y = y;

            if obj.rescale_and_log
                xb = [min(x), max(x)];
                obj.x_transform = @(x) rescale_and_log_x(x, xb);
            else
            	obj.x_transform = @(x) x;
            end

			[y_hat, iu] = unique(obj.keval(obj.x));
			obj.inv_interp = griddedInterpolant(y_hat,...
				obj.x(iu), 'spline', 'nearest');
		end

		% function set_reflection(obj)
		% 	preflect0 = 
		% end

		% function fit_polys(obj)
		% 	xlower = (obj.y <= 0.1);
		% 	obj.p_lower = polyfit(obj.x(xlower), obj.y(xlower), 3);
		% end

		function y_out = keval(obj, x_query)
			if numel(obj.force_fit_cdf_low) == 2
				adj_range = (obj.y >= obj.force_fit_cdf_low(1)) ...
					& (obj.y <= obj.force_fit_cdf_low(2));

				adjustment = ones(size(obj.x));
				adjustment(adj_range) = linspace(0.1, 1, sum(adj_range));
				h_adj = adjustment .* obj.h;
			else
				h_adj = obj.h;
			end
			
			y_out = zeros(size(x_query));
			x_qtrans = obj.x_transform(reshape(x_query, 1, []));
			x_trans = obj.x_transform(obj.x);

			ii = 1;
			for xq = x_qtrans
				v = abs(xq - x_trans) ./ h_adj;

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

				kd(x_trans - xq > xq) = 0;

				if sum(kd) > 0
					y_out(ii) = sum(obj.y .* kd) ./ sum(kd);
				else
					y_out(ii) = obj.y(1);
				end
				ii = ii + 1;
			end
        end
        
        function x_hat = keval_inv(obj, y_query)
        	x_hat = obj.inv_interp(y_query);
        end

        function ax = make_plots(obj, plot_raw, plot_fitted)
        	if nargin == 1
        		plot_raw = true;
        		plot_fitted = true;
        	end

        	legends = {};

        	if plot_raw
	        	plot(obj.x, obj.y)
	        	legends = [legends 'Raw data'];
	        	hold on
	        end

	        if plot_fitted
	            plot(obj.x, obj.keval(obj.x));
	            legends = [legends 'Kernel-smoothed (Gaussian)'];
	        end

	        legend(legends{:})
            hold off

            ax = gca();
        end
	end
end

function x_scaled = rescale_and_log_x(x, xbounds)
	x_scaled = (x - xbounds(1)) / (xbounds(2) - xbounds(1));
	x_scaled = log(1 + 100 * x);
end