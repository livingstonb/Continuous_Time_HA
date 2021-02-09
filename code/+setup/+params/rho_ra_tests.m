function [outparams, n] = rho_ra_tests(param_opts)
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
    shared_params.income_dir = 'continuous_a/no_measurement_error';
    shared_params.r_b = 0.02 / 4;
    shared_params.transfer = 0;
    shared_params.Bequests = true;
    
    params = {};
    
    %% NO TRANSITORY INCOME RISK
    % Iterate over r_a, rho
    no_trans_inc_risk = shared_params;
    no_trans_inc_risk.no_transitory_incrisk = true;
    no_trans_inc_risk.income_dir = 'continuous_a/no_measurement_error';
    
    no_trans_inc_risk.kappa1 = 5;
    no_trans_inc_risk.kappa2 = 2;
    
    rhos = 0.0009:0.0002:0.0022;
    r_as = 0.00505:0.00005:0.006;
    
    ii = 1;
    for rho = rhos
        for r_a = r_as
            params = [params no_trans_inc_risk];
            params{ii}.rho = rho;
            params{ii}.r_a = r_a;
            ii = ii + 1;
        end
    end
    
    %% OTHER
    % Iterate over r_a, rho
    high_kappa1_mid_kappa2 = shared_params;
    high_kappa1_mid_kappa2.no_transitory_incrisk = false;
    high_kappa1_mid_kappa2.income_dir = 'continuous_a/measurement_error_20pc';
    
    high_kappa1_mid_kappa2.kappa1 = 5;
    high_kappa1_mid_kappa2.kappa2 = 0.5;
    
    rhos = 0.0009:0.0002:0.003;
    r_as = 0.00505:0.002:0.03;
    
    for rho = rhos
        for r_a = r_as
            params = [params high_kappa1_mid_kappa2];
            params{ii}.rho = rho;
            params{ii}.r_a = r_a;
            ii = ii + 1;
        end
    end
   
    
    %% DO NOT CHANGE THIS SECTION
    n = numel(params);
    outparams = params{param_opts.param_index};
end