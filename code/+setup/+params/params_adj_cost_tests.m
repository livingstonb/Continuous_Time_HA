function outparams = params_adj_cost_tests(runopts)
	import HACTLib.aux.set_shared_fields

	chi1s = [0.01, 0.05, 0.1, 0.15, 0.2];
    chi2s = [0.1, 0.25, 0.5];

    ii = 1;

    shared_params.nb = 40;
    shared_params.nb_KFE = 40;
    shared_params.na = 40;
    shared_params.na_KFE = 40;
    % shared_params.min_grid_spacing = -1e10;
    shared_params.bmax = 20;
    shared_params.amax = 50;
    shared_params.b_gcurv_pos = 0.3;
    shared_params.a_gcurv = 0.3;
    shared_params.OneAsset = 0;
    shared_params.income_dir = 'continuous_b';
    shared_params.chi0 = 0;
    shared_params.a_lb = 0.25;

    %% --------------------------------------------------------------------
    % BASELINE
    % ---------------------------------------------------------------------
    params(ii).name = 'baseline';
    params(ii).chi1 = 0.15;
    params(ii).chi2 = 0.25;
    params(ii).rho = 0.015440584992491;
    params(ii).r_a = 0.0190643216;
    ii = ii + 1;

    %% --------------------------------------------------------------------
    % TARGET MEDIAN WEALTH STATS
    % ---------------------------------------------------------------------
    for chi1 = chi1s
    	for chi2 = chi2s
		    params(ii).name = sprintf('chi1=%g, chi2=%g',chi1, chi2);
		    params(ii).chi1 = chi1;
		    params(ii).chi2 = chi2;
		    params(ii).rho = 0.015440584992491;
		    params(ii).r_a = 0.02;
		    ii = ii + 1;
		end
	end

	params = set_shared_fields(params, shared_params);

    %% DO NOT CHANGE THIS SECTION
    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts,chosen_param);

    %% ATTACH CALIBRATOR
    [fn_handle, x0] = median_wealth_calibrator(outparams, runopts);
    outparams.set("calibrator", fn_handle, true);
    outparams.set("x0_calibration", x0, true);
end

function [fn_handle, x0] = median_wealth_calibrator(p, runopts)
	import HACTLib.model_objects.AltCalibrator

	% Vary (rho, r_a) to match median(a+b) = 1.6, median(b) = 0.1
	param_name = {'rho', 'r_a'};
	stat_name = {'median_totw', 'median_liqw'};
	stat_target = [1.6, 0.1];
	inits = [p.rho, p.r_a];

	rho_bounds = [0.002, 0.4];
	ra_bounds = [p.r_b+1e-4, 0.3];

	calibrator = AltCalibrator(p, runopts, param_name,...
        stat_name, stat_target);

    calibrator.set_param_bounds(1, rho_bounds);
	calibrator.set_param_bounds(2, ra_bounds);

    for ii = 1:numel(inits)
    	x0(ii) = calibrator.convert_to_solver_input(inits(ii), ii);
    end

    fn_handle = @(x) calibrator.fn_handle(x, p);
end