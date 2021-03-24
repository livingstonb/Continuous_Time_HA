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
    
    incomedirs = {'continuous_a/no_measurement_error',...
        'continuous_a/measurement_error_20pc',...
        'continuous_a/measurement_error_33pc',...
        'continuous_a/measurement_error_50pc'};

    IncomeDescriptions = {'cont_a, no meas err',...
        'cont_a, meas err 20pc',...
        'cont_a, meas err 33pc',...
        'cont_a, meas err 50pc'};

    experiment = true;
    if experiment
        iy = 1;
        params = shared_params;
        params.calibration_vars = {'rho', 'r_b'};
        params.calibration_stats = {'median_totw', 'median_liqw'};
        params.calibration_targets = [1.54, 0.05];
        params.calibration_scales = [1, 100];
        params.income_dir = incomedirs{iy};
        params.IncomeDescr = IncomeDescriptions{iy};
        param_opts.param_index = 1;

        params.kappa1 = 1;
        params.kappa2 = 1;
        params.rho = 0.009;
        params.r_b = 0.005;
        params.r_a = 0.009;

        rho_bds = [0.007, 0.015];
        r_b_bds = [0.004, 0.008];
        params.KFE_maxiters = 3e5;

        % Set calibrator
        params.calibration_bounds = {rho_bds, r_b_bds};
        params.calibration_backup_x0 = {};

        params = {params};
    else
        params = {};

        %% TARGET MEDIAN TOTAL WEALTH AND MEDIAN LIQUID WEALTH
        % Iterate over r_a, rho
        median_calibration = shared_params;
        median_calibration.calibration_vars = {'rho', 'r_a'};

        kappa_1s = [0.2:0.2:1, 1.5:0.5:5];
        kappa_2s = [0.5, 1.0, 1.5];

        calibrations = {median_calibration};

        ii = 1;
        group_num = 0;
        for icalibration = [1]
            for kappa2 = kappa_2s
                for iy = 1:3
                    group_num = group_num + 1;
                    for kappa1 = kappa_1s
                        params = [params {calibrations{icalibration}}];
                        params{ii}.name = sprintf('iy=%d, kappa2=%g', iy, kappa2);
                        params{ii}.kappa2 = kappa2;
                        params{ii}.income_dir = incomedirs{iy};
                        params{ii}.IncomeDescr = IncomeDescriptions{iy};
                        params{ii}.group_num = group_num;
                        
                        if params{ii}.no_transitory_incrisk
                            % params{ii}.rho = 0.001;
                            % params{ii}.r_a = 0.0052;
                            % params{ii}.calibration_bounds = {[0.0008, 0.003],...
                            % [shared_params.r_b + 0.0003, 0.009]};
                            % params{ii}.calibration_backup_x0 = {};
                        else
                            params{ii}.kappa1 = kappa1;
                            params{ii}.rho = 0.01;
                            params{ii}.r_a = 0.015;

                            rho_bds = [0.0075, 0.018];
                            r_a_bds = [0.006, 0.023];
                            params{ii}.KFE_maxiters = 3e5;

                            % params{ii}.rho = mean(rho_bds);
                            % params{ii}.r_a = mean(r_a_bds);

                            % Set calibrator
                            params{ii}.calibration_bounds = {rho_bds, r_a_bds};
                            params{ii}.calibration_backup_x0 = {};
                        end
                        params{ii}.calibration_stats = {'median_totw', 'median_liqw'};
                        params{ii}.calibration_targets = [1.54, 0.05];
                        params{ii}.calibration_scales = [1, 100];

                        ii = ii + 1;
                    end
                end
            end
        end
    end
    
    %% DO NOT CHANGE THIS SECTION
    n = numel(params);
    outparams = params{param_opts.param_index};
end