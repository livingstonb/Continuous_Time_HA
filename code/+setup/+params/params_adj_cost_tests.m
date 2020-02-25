function outparams = params_adj_cost_tests(runopts)
	import HACTLib.aux.set_shared_fields

	% chi1s = [0.01, 0.05, 0.1, 0.15, 0.2];
 %    chi2s = [0.1, 0.25, 0.5];

    chi1s = [0.05 0.1 0.2:0.2:1];
    chi2s = [0.01:0.01:0.05 0.1:0.1:0.5];
    ii = 1;

    shared_params.nb = 40;
    shared_params.nb_KFE = 40;
    shared_params.na = 40;
    shared_params.na_KFE = 40;
    shared_params.min_grid_spacing = -1e10;
    shared_params.mpc_shocks_dollars = [-1, -500, -5000, 1, 500, 5000];
    shared_params.mpc_shocks = shared_params.mpc_shocks_dollars ./ 72000;
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
    % Secondary targets are total HtM 30-35% and PHtM ~10%
    % i.e. about 1/3 is HtM and 2/3 of those are wealthy
    %
    % Need P(b <= 1/6 quarterly inc) ~= 1/3
    % Need Wealthy HtM / Total HtM ~= 2/3

    for chi1 = chi1s
    	for chi2 = chi2s
		    params(ii).name = sprintf('chi1=%g, chi2=%g',chi1, chi2);
		    params(ii).chi1 = chi1;
		    params(ii).chi2 = chi2;
		    params(ii).rho = 0.01;
		    params(ii).r_a = 0.015;
            params(ii).r_b = 0;
		    ii = ii + 1;
		end
	end

	params = set_shared_fields(params, shared_params);

    %% DO NOT CHANGE THIS SECTION
    % Use runopts.param_index to choose which specification to select
    chosen_param = params(runopts.param_index);

    % Create Params object
    outparams = HACTLib.model_objects.Params(runopts, chosen_param);

    %% ATTACH CALIBRATOR
    if runopts.calibrate
        if strcmp(outparams.name, 'baseline')
        	[fn_handle, x0] = mean_wealth_calibrator(outparams, runopts);
        else
    	    [fn_handle, x0] = median_wealth_calibrator(outparams, runopts);
    	end

    	outparams.set("calibrator", fn_handle, true);
        outparams.set("x0_calibration", x0, true);
    end
end

function [fn_handle, x0] = mean_wealth_calibrator(p, runopts)
	import HACTLib.model_objects.AltCalibrator

	param_name = {'rho', 'r_a'};
	stat_name = {'totw', 'liqw'};
	stat_target = [3.5, 0.5];
	inits = [p.rho, p.r_a];

	rho_bounds = [0.002, 0.2];
	ra_bounds = [p.r_b+1e-3, 0.1];

	calibrator = AltCalibrator(p, runopts, param_name,...
        stat_name, stat_target);

    calibrator.set_param_bounds(1, rho_bounds);
	calibrator.set_param_bounds(2, ra_bounds);

    for ii = 1:numel(inits)
    	x0(ii) = calibrator.convert_to_solver_input(inits(ii), ii);
    end

    fn_handle = @(x) calibrator.fn_handle(x, p);
end

function [fn_handle, x0] = median_wealth_calibrator(p, runopts)
	import HACTLib.model_objects.AltCalibrator

	% Vary (rho, r_a) to match median(a+b) = 1.6, median(b) = 0.1
	param_name = {'rho', 'r_a'};
	stat_name = {'median_totw', 'median_liqw'};
	stat_target = [1.7, 0.1];
	inits = [p.rho, p.r_a];

	rho_bounds = [0.002, 0.2];
    ra_bounds = [p.r_b+1e-3, 0.1];

	calibrator = AltCalibrator(p, runopts, param_name,...
        stat_name, stat_target);

    calibrator.set_param_bounds(1, rho_bounds);
	calibrator.set_param_bounds(2, ra_bounds);

    for ii = 1:numel(inits)
    	x0(ii) = calibrator.convert_to_solver_input(inits(ii), ii);
    end

    fn_handle = @(x) calibrator.fn_handle(x, p);
end

% function [fn_handle, x0] = four_stat_calibrator(p, runopts)
%     import HACTLib.model_objects.AltCalibrator

%     % Vary (rho, r_a) to match median(a+b) = 1.6, median(b) = 0.1
%     % and (chi1, chi2) to match HtM stats
%     param_name = {'rho', 'r_a', 'chi1', 'chi2'};
%     stat_name = {'median_totw', 'median_liqw',...
%         'ratio_WHtM_HtM_sixth', 'HtM_one_sixth_Q_lwealth'};
%     stat_target = [1.7, 0.1, 0.66, 0.33];
%     inits = [p.rho, p.r_a, 0.15, 0.25];

%     rho_bounds = [0.001, 0.1];
%     ra_bounds = [p.r_b+1e-4, 0.3];
%     chi1_bounds = [0.01, 5];
%     chi2_bounds = [0.01, 5];

%     calibrator = AltCalibrator(p, runopts, param_name,...
%         stat_name, stat_target);

%     calibrator.set_param_bounds(1, rho_bounds);
%     calibrator.set_param_bounds(2, ra_bounds);
%     calibrator.set_param_bounds(3, chi1_bounds);
%     calibrator.set_param_bounds(4, chi2_bounds);

%     for ii = 1:numel(inits)
%         x0(ii) = calibrator.convert_to_solver_input(inits(ii), ii);
%     end

%     fn_handle = @(x) calibrator.fn_handle(x, p);
% end