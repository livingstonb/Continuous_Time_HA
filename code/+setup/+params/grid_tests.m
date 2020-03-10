function outparams = table_tests(param_opts, param_index)
    % Create structure array 'params', and output a Params instance
    % of the structure in the 'index' entry, i.e. 1,2,3,.

    import HACTLib.aux.set_shared_fields

    shared_params = param_opts;
    shared_params.chi0 = 0;
    shared_params.chi1 = 0.15;
    shared_params.chi2 = 0.25;
    shared_params.a_lb = 0.25;
    shared_params.OneAsset = 0;
    shared_params.income_dir = 'continuous_b';
    shared_params.rho = 0.015;
    shared_params.r_a = 0.02;
    shared_params.calibration_vars = {'rho', 'r_a'};
    shared_params.calibration_stats = {'totw', 'liqw'};
    shared_params.calibration_targets = [3.5, 0.5];
    shared_params.calibration_bounds = {[0.001, 0.03], [0.008, 0.04]};
    shared_params.min_grid_spacing = 1e-5;
    shared_params.saveGrids = true;

    nxs = 50;
    curvs = [0.2, 0.3];
    maxes = {[5, 10], [25, 50], [50, 100]};

    ii = 1;
    for nx = nxs
        for curv = curvs
            for imax = 1:3
                params{ii} = shared_params;
                params{ii}.name = sprintf('test %d', ii);
                params{ii} = set_grid(params{ii}, nx, curv)
                params{ii} = set_maxes(params{ii}, maxes{imax});
                ii = ii + 1;
            end
        end
    end

    maxes = [10, 30];
    params{ii} = shared_params;
    params{ii}.name = 'brute force, lower maxes';
    params{ii} = set_grid(params{ii}, 200, 0.3)
    params{ii} = set_maxes(params{ii}, maxes);
    ii = ii + 1;

    maxes = [20, 60];
    params{ii} = shared_params;
    params{ii}.name = 'brute force, higher maxes';
    params{ii} = set_grid(params{ii}, 220, 0.3)
    params{ii} = set_maxes(params{ii}, maxes);

    % ---------------------------------------------------------------------
    % Use param_index to choose which specification to select
    outparams = params{param_index};
end

function params_in = set_grid(params_in, nx, curv)
    params_in.nb = nx;
    params_in.nb_KFE = nx;
    params_in.na = nx;
    params_in.na_KFE = nx;
    params_in.a_gcurv = curv;
    params_in.b_gcurv_pos = curv;
end

function params_in = set_maxes(params_in, maxes)
    params_in.bmax = maxes(1);
    params_in.amax = maxes(2);
end