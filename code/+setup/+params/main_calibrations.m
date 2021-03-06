function [outparams, n] = main_calibrations(param_opts)
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
    
    params = {};
    
    %% TARGET MEDIAN TOTAL WEALTH AND MEDIAN LIQUID WEALTH
    % Iterate over r_a, rho
    median_calibration = shared_params;
    median_calibration.calibration_vars = {'rho', 'r_a'};
    
    % kappa_0s = [0];
    % kappa_1s = [0.025:0.025:0.1 0.25:0.25:1.75 2:1:10];
    % kappa_2s = [0.1, 0.5, 0.75, 1.0, 2.0];

    kappa_0s = [0];
    kappa_1s = [1];
    kappa_2s = [1];

    anninc = shared_params.numeraire_in_dollars;
    a_lbs_dollars = [100, 500, 1000, 2000, 5000, 10000, 20000];
    a_lbs = a_lbs_dollars / anninc;
    % a_lbs = [a_lbs, 0.25];

    % Since ra << kappa1 ^ (-1 / kappa2) to prevent illiquid
    % asset overaccumulation at the top, want to avoid
    % case of high kappa1, low kappa 2. Split each param
    % into two groups: low and high
    % kappa1s = {[0.5, 1], [2, 5, 10]};
    % kappa2s = {[0.1, 0.25, 0.5], [1, 1.25]};
    % icombinations = {[1, 1], [1, 2], [2, 2]};
    % kappa1s = {1:0.25:3};
    % kappa2s = {[0.025 0.05 0.075 0.1:0.1:0.5]};
    % icombinations = {[1, 1]};
    
    calibrations = {median_calibration};
    calibration_labels = {'MEDIAN TARGETS'};

%     incomedirs = {'results_feb_2020_3pt/continuous_a/no_measurement_error',...
%         'results_feb_2020_3pt/continuous_a/measurement_error_50pc'};
% 
%     IncomeDescriptions = {'cont_a 2/20, no meas err',...
%         'cont_a 2/20, meas err 50pc'};

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
        for iy = 2
            for a_lb_dollar = a_lbs_dollars
                a_lb = a_lb_dollar / anninc;
                for kappa0 = kappa_0s
                    for kappa1 = kappa_1s
                        for kappa2 = kappa_2s
                            params = [params {calibrations{icalibration}}];
                            % params{ii}.name = sprintf('%s, kappa0=%g, kappa1=%g, kappa2=%g',...
                            %     calibration_labels{icalibration}, kappa0, kappa1, kappa2);
                            params{ii}.name = sprintf('a_lb=%g, kappa1=%g, kappa2=%g',...
                                a_lb, kappa1, kappa2);
                            params{ii}.kappa0 = kappa0;
                            params{ii}.kappa1 = kappa1;
                            params{ii}.kappa2 = kappa2;
                            params{ii}.a_lb = a_lb;
                            params{ii}.income_dir = incomedirs{iy};
                            params{ii}.IncomeDescr = IncomeDescriptions{iy};
                            
                            if params{ii}.no_transitory_incrisk
                                params{ii}.rho = 0.001;
                                params{ii}.r_a = 0.0052;
                                params{ii}.calibration_bounds = {[0.0008, 0.003],...
                                [shared_params.r_b + 0.0003, 0.009]};
                                params{ii}.calibration_backup_x0 = {};
                            else
                                printbds = (ii == param_opts.param_index);
                                % [rho_bds, r_a_bds] = setup.params.get_rho_ra_bounds(kappa1, kappa2, printbds);
                                [rho_bds, r_a_bds] = setup.params.get_rho_ra_bounds_a_lb_runs(...
                                    a_lb_dollar, printbds);
                                params{ii}.KFE_maxiters = 3e5;
        
                                params{ii}.rho = mean(rho_bds);
                                params{ii}.r_a = mean(r_a_bds);

                                % Set calibrator
                                params{ii}.calibration_bounds = {rho_bds, r_a_bds};
                                params{ii}.calibration_backup_x0 = {};
                            end
                            params{ii}.calibration_stats = {'median_totw', 'median_liqw'};
                            params{ii}.calibration_targets = [scf.median_totw, scf.median_liqw];
                            params{ii}.calibration_scales = [1, 100];
        
                            ii = ii + 1;
                        end
                    end
                end
            end
        end
    end
    
    %% DO NOT CHANGE THIS SECTION
    n = numel(params);
    outparams = params{param_opts.param_index};
end