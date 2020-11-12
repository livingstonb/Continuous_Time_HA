classdef HACTCalibrator < EconToolsML.Calibrator
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	methods
		function obj = HACTCalibrator(params, variables,...
			target_names, target_values)

			obj = obj@EconToolsML.Calibrator(params, variables, target_names, target_values);
			
			obj.main_handle = @(curr_params) main(curr_params, 'iterating', true);
		end

		function construct_options_struct(obj, params)
			obj.options = struct();
			obj.options.ComputeMPCS = params.ComputeMPCS;
			obj.options.ComputeMPCS_illiquid = params.ComputeMPCS_illiquid;
			obj.options.SimulateMPCS = params.SimulateMPCS;
			obj.options.ComputeMPCS_news = params.ComputeMPCS_news;
			obj.options.SimulateMPCS_news = params.SimulateMPCS_news;
			obj.options.NoRisk = params.NoRisk;
		end

		function value = get_results_value(obj, results, variable_name)
			value = results.(variable_name).value;
		end

		function dv = adjust_dv(obj, results, current_params, dv)
			if ismember('r_a', obj.target_names)
                % If illiquid wealth is very close to zero, introduce an
                % ad-hoc reward on increasing the illiquid return. The
                % effect of this vanishes as illiquid wealth gets larger.
				z = results.stats.illiqw.value - 0.01;

				c = min(max(z, 0), 0.5);
				m = (1 + cos(c * pi * 2)) / 10;
				dv(end+1) = (0.01 / (0.2 + current_params.r_a)) * m;
			else
				% Do nothing
				dv = dv;
			end
		end
	end
end