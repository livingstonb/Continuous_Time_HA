classdef KernelSmoother < handle
	properties
		x;
		inv_interp;
		y;
		h;
		ktype = 'gaussian';

		force_fit_bottom = true;

		log_transform_const = 0.1;
		log_transform = false;
		transform_x = @(x) x;
	end

	methods
		function obj = KernelSmoother(ktype, log_transform,...
			log_transform_const)
			if nargin == 1
				obj.ktype = ktype;
			end

			if nargin >= 2
				obj.ktype = ktype;
				obj.log_transform = log_transform;

				if nargin == 3
					obj.log_transform_const = log_transform_const;
				end
			end

			if obj.log_transform
				obj.transform_x = ...
					@(x) exp(x) - obj.log_transform_const;
			end
		end

		function set_for_cdf_smoothing(obj, pmf, values, h)
			arr_sorted = sortrows([values(:), pmf(:)]);

			[v_u, iu] = unique(arr_sorted(:,1), 'last');
			cdf_u = cumsum(arr_sorted(:,2));
			cdf_u = cdf_u(iu);

			if obj.log_transform
				v_u = log(obj.log_transform_const + v_u);
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
				obj.x = log(obj.log_transform_const + x);
			else
				obj.x = x;
			end
			obj.y = y;
			obj.h = h;

			[y_hat, iu] = impose_monotonic(obj.transform_x(obj.x),...
				obj.keval(obj.transform_x(obj.x)));

			obj.inv_interp = griddedInterpolant(y_hat,...
				obj.transform_x(obj.x(iu)),...
        		'spline', 'nearest');
		end

		function x_out = keval(obj, x_query)
			if obj.log_transform
				x_query = log(obj.log_transform_const + x_query);
			end

			if obj.force_fit_bottom
				adj_range = (obj.y >= 0) & (obj.y <= 0.05);

				adjustment = ones(size(obj.x));
				adjustment(adj_range) = linspace(0.01, 1, sum(adj_range))';
				h_adj = adjustment .* obj.h;
			else
				h_adj = obj.h;
			end

			x_query = reshape(x_query, 1, []);
			v = abs(x_query - obj.x) ./ h_adj;

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
        	x_hat = obj.inv_interp(y_query);
        end

        function make_plots(obj, plot_raw, plot_fitted)
        	x_grid = obj.transform_x(obj.x);

        	if nargin == 1
        		plot_raw = true;
        		plot_fitted = true;
        	end

        	legends = {};

        	if plot_raw
	        	plot(x_grid, obj.y)
	        	legends = [legends 'Raw data'];
	        	hold on
	        end

	        if plot_fitted
	            plot(x_grid, obj.keval(obj.transform_x(obj.x)));
	            legends = [legends 'Kernel-smoothed (Gaussian)'];
	        end

	        legend(legends{:})
            hold off
        end
        
        function plot_cdf(obj, plot_raw, plot_fitted)
        	if nargin == 1
        		plot_raw = true;
        		plot_fitted = true;
        	end

        	if plot_raw
            	plot(obj.transform_x(obj.x), cumsum(obj.y))
            	legend('Raw data')
            end

            if plot_fitted
            	if plot_raw
            		hold on
            	end
	            plot(obj.transform_x(obj.x), cumsum(obj.keval(obj.x)));
	            legend('Raw data', 'Kernel-smoothed (Gaussian)')
	            hold off
           end
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

function [arr_out, iu] = impose_monotonic(agrid, arr_in)
	increasing = false(size(arr_in));

	last_val = arr_in(1);
	end_val = arr_in(end);
	increasing(1) = true;
	increasing(end) = true;
	for ia = 2:numel(arr_in)-1
		if (arr_in(ia) > last_val) && (arr_in(ia) < end_val)
			increasing(ia) = true;
			last_val = arr_in(ia);
		end
	end

	arr_out = arr_in(increasing);
	iu = find(increasing);
end