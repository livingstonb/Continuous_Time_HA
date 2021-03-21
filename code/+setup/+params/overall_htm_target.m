function [outparams, n] = overall_htm_target(param_opts)
	import setup.params.scf2019struct
    
    scf = scf2019struct();

	shocks = [-1, -500, -5000, 1, 500, 5000];

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / (scf.quarterly_earnings * 4);
    shared_params.numeraire_in_dollars = (scf.quarterly_earnings * 4);
    shared_params.nb = 50;
    shared_params.nb_KFE = 50;
    shared_params.na = 50;
    shared_params.na_KFE = 50;
    
    shared_params.bgrid_term1_weight = 0.01;
    shared_params.bgrid_term1_curv = 0.8;
    shared_params.agrid_term1_weight = 0.01;
    shared_params.agrid_term1_curv = 0.6;

    shared_params.b_gcurv_pos = 0.1;
    shared_params.a_gcurv = 0.2;

    shared_params.a_lb = 0.25;
    shared_params.bmax = 25;
    shared_params.amax = 200;
    shared_params.OneAsset = 0;
    shared_params.income_dir = 'continuous_a';
    shared_params.r_b = 0.02 / 4;
    shared_params.transfer = 0;
    shared_params.Bequests = true;
    shared_params.no_transitory_incrisk = false;

    shared_params.kappa0 = 0;

    anninc = shared_params.numeraire_in_dollars;
    shared_params.a_lb = 500 / anninc;
    
    params = {};
    
    %% TARGET MEDIAN TOTAL WEALTH AND MEDIAN LIQUID WEALTH
    % Iterate over r_a, rho
    median_calibration = shared_params;
    median_calibration.calibration_vars = {'rho', 'r_a', 'kappa1'};

    kappa_2s = [0.5, 1.0, 1.5];

    calibrations = {median_calibration};
    calibration_labels = {'MEDIAN TARGETS'};

    incomedirs = {'continuous_a/no_measurement_error',...
        'continuous_a/measurement_error_20pc',...
        'continuous_a/measurement_error_33pc',...
        'continuous_a/measurement_error_50pc'};

    IncomeDescriptions = {'cont_a, no meas err',...
        'cont_a, meas err 20pc',...
        'cont_a, meas err 33pc',...
        'cont_a, meas err 50pc'};
    
    ii = 1;
    for icalibration = [1]
        for kappa2 = kappa_2s
            for iy = 1:3
                params = [params {calibrations{icalibration}}];
                params{ii}.name = sprintf('iy=%d, kappa2=%g', iy, kappa2);
                params{ii}.kappa2 = kappa2;
                params{ii}.income_dir = incomedirs{iy};
                params{ii}.IncomeDescr = IncomeDescriptions{iy};
                
                if params{ii}.no_transitory_incrisk
                    params{ii}.rho = 0.001;
                    params{ii}.r_a = 0.0052;
                    params{ii}.calibration_bounds = {[0.0008, 0.003],...
                    [shared_params.r_b + 0.0003, 0.009]};
                    params{ii}.calibration_backup_x0 = {};
                else
                    params{ii}.kappa1 = 0.5;
                    params{ii}.rho = 0.00875;
                    params{ii}.r_a = 0.015;

                    kappa1_bds = [0.2, 5];
                    rho_bds = [0.005, 0.015];
                    r_a_bds = [0.007, 0.02];
                    params{ii}.KFE_maxiters = 3e5;

                    % params{ii}.rho = mean(rho_bds);
                    % params{ii}.r_a = mean(r_a_bds);

                    % Set calibrator
                    params{ii}.calibration_bounds = {rho_bds, r_a_bds, kappa1_bds};
                    params{ii}.calibration_backup_x0 = {};
                end
                params{ii}.calibration_stats = {'median_totw', 'median_liqw', 'liqw_lt_ysixth'};
                params{ii}.calibration_targets = [1.54, 0.5, 0.39];
                params{ii}.calibration_scales = [1, 100, 5];

                ii = ii + 1;
            end
        end
    end
    
    %% DO NOT CHANGE THIS SECTION
    n = numel(params);
    outparams = params{param_opts.param_index};
end