function outparams = params_adj_cost_tests(param_opts, param_index)
	import HACTLib.aux.set_shared_fields

	shocks = [-1, -500, -5000, 1, 500, 5000];
    dollars = 72000;

    shared_params = param_opts;
    shared_params.mpc_shocks = shocks / dollars;
    shared_params.numeraire_in_dollars = dollars;
    shared_params.nb = 50;
    shared_params.nb_KFE = 50;
    shared_params.na = 50;
    shared_params.na_KFE = 50;
    shared_params.min_grid_spacing = 0;
    shared_params.b_gcurv_pos = 0.2;
    shared_params.a_gcurv = 0.2;
    shared_params.a_glinear = 0.02;
    shared_params.b_glinear = 0.01;
    shared_params.a_lb = 0.25;
    shared_params.bmax = 20;
    shared_params.amax = 75;
    shared_params.OneAsset = 0;
    shared_params.income_dir = 'continuous_b';
    % shared_params.calibration_vars = {'rho', 'r_a'};
    % shared_params.calibration_bounds = {[0.001, 0.05], [0.006, 0.1]};

    shared_params.calibration_vars = {'rho', 'r_a'};
    shared_params.calibration_bounds = {[0.001, 0.05], [0.006, 0.1]};
    shared_params.calibration_stats = {'totw', 'liqw'};
    shared_params.calibration_targets = [3.5, 0.5];

    ii = 1;

    %% --------------------------------------------------------------------
    % BASELINE
    % ---------------------------------------------------------------------
    % params{ii} = shared_params;
    % params{ii}.name = 'baseline';
    % params{ii}.chi0 = 0;
    % params{ii}.a_lb = 0.25;
    % params{ii}.chi1 = 0.15;
    % params{ii}.chi2 = 0.25;
    % params{ii}.rho = 0.0144188;
    % params{ii}.r_a = 0.0188616;
    % params{ii}.calibration_stats = {'totw', 'liqw'};
    % params{ii}.calibration_targets = [3.5, 0.5];
    % ii = ii + 1;

    %% --------------------------------------------------------------------
    % TARGET MEDIAN WEALTH STATS
    % ---------------------------------------------------------------------
    % Secondary targets are total HtM 30-35% and PHtM ~10%
    % i.e. about 1/3 is HtM and 2/3 of those are wealthy
    %
    % Need P(b <= 1/6 quarterly inc) ~= 1/3
    % Need Wealthy HtM / Total HtM ~= 2/3

    % chi1s = [0.1 0.5 1];
    % chi2s = [0.05 0.1 0.5 1];
    % chi0s = [0 0.01 0.1];

    rho0s = [0.01, 0.01];
    ra_0s = [0.015, 0.02];

    kappa_0s = [0, 0.03];
    kappa_1s = [1, 2, 5];
    kappa_2s = [0.5, 1, 1.5];

    for rho = rho0s
        for r_a = ra_0s
            for kappa0 = kappa_0s
                for kappa1 = kappa_1s
                    for kappa2 = kappa_2s
                        params{ii} = shared_params;
                        params{ii}.name = sprintf('kappa0=%g, kappa1=%g, kappa2=%g',...
                            kappa0, kappa1, kappa2);
                        params{ii}.kappa0 = kappa0;
                        params{ii}.kappa1 = kappa1;
                        params{ii}.kappa2 = kappa2;
                        params{ii}.rho = rho;
                        params{ii}.r_a = r_a;
                        params{ii}.r_b = 0.005;

                        % params{ii}.KFE_maxiters = 1e5;
                        % params{ii}.calibration_vars = {'rho', 'r_a'};
                        % params{ii}.calibration_bounds = {[0.001, 0.02], [0.006, 0.1]};
                        % params{ii}.calibration_stats = {'median_totw', 'median_liqw'};
                        % params{ii}.calibration_targets = [1.7, 0.1];
                        ii = ii + 1;
                    end
                end
            end
        end
    end

    %% DO NOT CHANGE THIS SECTION
    % Use param_index to choose which specification to select
    outparams = params{param_index};
end