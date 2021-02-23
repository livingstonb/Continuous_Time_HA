function [outparams, n] = kappa2_075_runs(param_opts)
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
    shared_params.income_dir = 'continuous_a/measurement_error_33pc';
    shared_params.r_b = 0.02 / 4;
    shared_params.transfer = 0;
    shared_params.Bequests = true;

    shared_params.kappa2 = 0.75;

    shared_params.calibration_stats = {'median_totw', 'median_liqw'};
    shared_params.calibration_targets = [scf.median_totw, scf.median_liqw];
    shared_params.calibration_scales = [1, 100];
    shared_params.calibration_backup_x0 = {};
    shared_params.calibration_vars = {'rho', 'r_a'};
    
    params = {};
    
    %% NO TRANSITORY INCOME RISK
    
    rhos = 0.0009:0.0002:0.0022;
    r_as = 0.00505:0.00005:0.006;

    incomedirs = {'continuous_a/no_measurement_error',...
        'continuous_a/measurement_error_20pc',...
        'continuous_a/measurement_error_33pc',...
        'continuous_a/measurement_error_50pc'};

    IncomeDescriptions = {'cont_a, no meas err',...
        'cont_a, meas err 20pc',...
        'cont_a, meas err 33pc',...
        'cont_a, meas err 50pc'};

    kappa1s = [5, 7, 8, 9, 10];

    ii = 1;
    for iy = 1:5
        for kappa1 = kappa1s
            params = [params shared_params];

            params{ii}.kappa1 = kappa1;

            if (iy == 5)
                params{ii}.no_transitory_incrisk = true;
                params{ii}.income_dir = incomedirs{1};
                params{ii}.IncomeDescr = 'No trans risk';
            else
                params{ii}.no_transitory_incrisk = false;
                params{ii}.income_dir = incomedirs{iy};
                params{ii}.IncomeDescr = IncomeDescriptions{iy};
            end

            if (kappa1 == 5)
                rho_bds = [0.011, 0.013];
                r_a_bds = [0.014, 0.015];
                params{ii}.rho = 0.012;
                params{ii}.r_a = 0.0145;
            elseif (kappa1 == 7)
                rho_bds = [0.011, 0.013];
                r_a_bds = [0.0133, 0.0148];
                params{ii}.rho = 0.0122;
                params{ii}.r_a = 0.0141;
            elseif (kappa1 == 8)
                rho_bds = [0.011, 0.013];
                r_a_bds = [0.013, 0.0145];
                params{ii}.rho = 0.012;
                params{ii}.r_a = 0.0137;
            elseif (kappa1 == 9)
                rho_bds = [0.011, 0.013];
                r_a_bds = [0.0125, 0.014];
                params{ii}.rho = 0.012;
                params{ii}.r_a = 0.0134;
            elseif (kappa1 == 10)
                rho_bds = [0.011, 0.013];
                r_a_bds = [0.012, 0.014];
                params{ii}.rho = 0.012;
                params{ii}.r_a = 0.013;
            end

            params{ii}.calibration_bounds = {rho_bds, r_a_bds};

            ii = ii + 1;
        end
    end
    
    %% DO NOT CHANGE THIS SECTION
    n = numel(params);
    outparams = params{param_opts.param_index};
end