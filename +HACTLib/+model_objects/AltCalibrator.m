classdef AltCalibrator < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties
		runopts;
		options;
		variables;
		target_names;
		target_values;
		nvars;

		lbounds;
		ubounds;

		x0;

		solver_handle;
	end

	methods
		function obj = AltCalibrator(params, runopts, variables,...
			target_names, target_values)
			obj.runopts = runopts;

			obj.options = struct();
			obj.options.ComputeMPCS = params.ComputeMPCS;
			obj.options.ComputeMPCS_illiquid = params.ComputeMPCS_illiquid;
			obj.options.SimulateMPCS = params.SimulateMPCS;
			obj.options.ComputeMPCS_news = params.ComputeMPCS_news;
			obj.options.SimulateMPCS_news = params.SimulateMPCS_news;
			obj.options.NoRisk = params.NoRisk;
			obj.options.SaveResults = params.SaveResults;

			obj.variables = variables;

			obj.target_names = target_names;
			obj.target_values = target_values;

			obj.nvars = numel(variables);
			obj.lbounds = NaN(obj.nvars, 1);
			obj.ubounds = NaN(obj.nvars, 1);

			for i_var = 1:obj.nvars
				obj.x0(i_var) = params.(obj.variables{i_var});
			end
		end

		function set_param_bounds(obj, i_var, bounds)
			obj.lbounds(i_var) = bounds(1);
			obj.ubounds(i_var) = bounds(2);
		end

		function set_handle(obj, p)
			obj.solver_handle = @(x) obj.fn_handle(x, p);
		end

		function dv = fn_handle(obj, x, current_params)
			obj.turn_off_param_options(current_params);

			for i_var = 1:obj.nvars
				current_params.set(obj.variables{i_var}, x(i_var));
			end
			stats = main(obj.runopts, current_params);

			for i_var = 1:obj.nvars
				
			end
			
			fprintf('\n\n---- For ')
			for i_var = 1:obj.nvars
				v(i_var) = stats.(obj.target_names{i_var});
				fprintf('%s = %g', obj.variables{i_var}, x(i_var))
				if i_var < obj.nvars
					fprintf(", ")
				else
					fprintf(':\n')
				end
			end
			dv = v - obj.target_values;

			fprintf('    Result was ')
			for i_var = 1:obj.nvars
				fprintf('%s = %g', obj.target_names{i_var}, v(i_var))
				if i_var < obj.nvars
					fprintf(", ")
				else
					fprintf('\n')
				end
			end

			obj.reset_param_options(current_params);
		end

		function turn_off_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, false, quiet);
			end
		end

		function reset_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, obj.options.(props{ip}), quiet);
			end
		end
	end
end
