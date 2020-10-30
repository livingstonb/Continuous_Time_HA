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
    shared_params.amax = 100;
    shared_params.rho = 0.012;
    shared_params.r_a = 0.015;
    shared_params.OneAsset = 0;
    shared_params.income_dir = 'continuous_b';
    
    params = {};
    
    %% TARGET MEDIAN TOTAL WEALTH AND MEDIAN LIQUID WEALTH
    % Iterate over r_a, rho
    median_calibration = shared_params;
    median_calibration.calibration_vars = {'rho', 'r_a'};
    median_calibration.calibration_bounds = {[0.001, 0.05], [0.0055, 0.05]};
    median_calibration.calibration_backup_x0 = {[0.004, 0.0065]};
    median_calibration.calibration_stats = {'median_totw', 'median_liqw'};
    median_calibration.calibration_targets = [scf.median_totw, scf.median_liqw];
    
    kappa_0s = [0, 0.01];

    % Since ra << kappa1 ^ (-1 / kappa2) to prevent illiquid
    % asset overaccumulation at the top, want to avoid
    % case of high kappa1, low kappa 2. Split each param
    % into two groups: low and high
    kappa1s = {[0.5, 1], [2, 5, 10]};
    kappa2s = {[0.1, 0.25, 0.5], [1, 1.25]};
    icombinations = {[1, 1], [1, 2], [2, 2]};
    
    calibrations = {median_calibration};
    calibration_labels = {'MEDIAN TARGETS'};
    
    ii = 1;
    for icalibration = [1]
    for kappa0 = kappa_0s
        for ic = 1:numel(icombinations)
            i_kappa1s = icombinations{ic}(1);
            i_kappa2s = icombinations{ic}(2);
            kappa1s_comb = kappa1s{i_kappa1s};
            kappa2s_comb = kappa2s{i_kappa2s};

            for kappa1 = kappa1s_comb
            for kappa2 = kappa2s_comb
                params = [params {calibrations{icalibration}}];
                params{ii}.name = sprintf('%s, kappa0=%g, kappa1=%g, kappa2=%g',...
                    calibration_labels{icalibration}, kappa0, kappa1, kappa2);
                params{ii}.kappa0 = kappa0;
                params{ii}.kappa1 = kappa1;
                params{ii}.kappa2 = kappa2;
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