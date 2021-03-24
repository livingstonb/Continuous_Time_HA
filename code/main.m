function stats = main(p, varargin)
    % Instantiates necessary classes and calls functions to solve the
    % model and compute statistics
    %
    % Parameters
    % ----------
    % p : a Params object containing model parameters
    %
    % Returns
    % -------
    % stats : a structure containing statistics from the solved model
    %
    % p : the Params object used to solve the model

    import HACTLib.model_objects.Grid
    import HACTLib.model_objects.Income

    % Parse options
    parser = inputParser;
    addOptional(parser, 'quiet', false);
    addOptional(parser, 'final', false);
    addOptional(parser, 'save_iteration_results', false);
    parse(parser, varargin{:});
    quiet = parser.Results.quiet;
    final = parser.Results.final;
    save_iteration_results = parser.Results.save_iteration_results;

    if quiet
        % Turn printing off
        fprintf_internal = @(varargin) fprintf('');
    else
        fprintf_internal = @(varargin) fprintf(varargin{:});
    end

    %% --------------------------------------------------------------------
    % CREATE GRID, INCOME OBJECTS
    % ---------------------------------------------------------------------
    % Main income process
    norisk = false;
	income = Income(fullfile('input', p.income_dir), p, norisk);

    % Turn off income risk (set y equal to the mean)
    norisk = true;
    income_norisk = Income('', p, norisk);

    % Update number of income points
    p.set("ny", income.ny, true);

    % Borrowing limit
    NBL = - min(income.y.vec + p.transfer) ...
            / (p.r_b_borr + p.deathrate * p.perfectannuities);
    if p.bmin <= -1e10
        % Set "loose" borrowing limit
        p.set("bmin", 0.95 * NBL, true);
    elseif p.bmin < 0
        % Check that borrowing limit does not violate NBL
        msg = sprintf('bmin < natural borrowing limit (%f)', NBL);
        assert(p.bmin > NBL, msg);
    end

	grd = Grid(p, income.ny, 'HJB').auto_construct(); % grid for HJB
    grd_norisk = Grid(p, 1, 'HJB').auto_construct();
	grdKFE = Grid(p, income.ny, 'KFE').auto_construct();% grid for KFE
    grdKFE_norisk = Grid(p, 1, 'KFE').auto_construct();

    % Add net income variables
    income.set_net_income(p, grd, grdKFE);
    income_norisk.set_net_income(p, grd_norisk, grdKFE_norisk);

    % Instantiate and solve Model instance
    model = HACTLib.model_objects.Model(p, grd, grdKFE, income, 'quiet', quiet);
    model.initialize();
    [~, KFE, Au] = model.solve();

    if p.SolveNoRisk == 1
        % Solve model without income risk
        model_nr = HACTLib.model_objects.Model(...
            p, grd_norisk, grdKFE_norisk, income_norisk);
        model_nr.initialize();
        [~, KFE_nr, Au_nr] = model_nr.solve();
    end

    % if p.makePlots
    %     HACTLib.plots.make_wealth_histograms(model);
    % end

    %% ----------------------------------------------------------------
    % COMPUTE STATISTICS
    % -----------------------------------------------------------------
    fprintf_internal('\nComputing statistics\n')    
    stats = HACTLib.Statistics(p, income, grdKFE, KFE);
    stats.compute_statistics();

    %% ----------------------------------------------------------------
    % COMPUTE MPCs
    % -----------------------------------------------------------------
    import HACTLib.computation.MPCs

    % MPCs out of shock in current period
    mpc_finder = MPCs(p, income, grdKFE, p.mpc_options);
    mpc_finder_illiquid = MPCs(p, income, grdKFE, p.mpc_options_illiquid);
    
    if p.ComputeMPCS
    	fprintf_internal('\nComputing MPCs out of an immediate shock...\n')
        mpc_finder.solve(KFE, Au, stats.pmf);
    end

    if p.ComputeMPCS_illiquid
        fprintf_internal('\nComputing illiquid MPCs out of an immediate shock...\n')
        mpc_finder_illiquid.solve(KFE, Au, stats.pmf);
    end
    
    stats.add_mpcs(mpc_finder);
    stats.add_mpcs(mpc_finder_illiquid);

    % MPCs out of news
    shocks = [4, 5, 6];
    import HACTLib.computation.MPCsNews
    trans_dyn_solver = MPCsNews(p, income, grdKFE, shocks, p.mpcs_news_options);

    if p.ComputeMPCS_news
    	fprintf_internal('Computing MPCs out of news...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf,mpc_finder.cum_con_baseline);
    elseif p.SimulateMPCS_news
        fprintf_internal('Iterating backward to find policy functions for future shock...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf);
    end
    stats.add_mpcs_news(trans_dyn_solver);

    %% ----------------------------------------------------------------
    % SIMULATE MPCs
    % -----------------------------------------------------------------
    shocks = [4,5,6];
    shockperiod = 0;
    mpc_simulator = HACTLib.computation.MPCSimulator(...
    	p, income, grdKFE, KFE, shocks, shockperiod, p.mpcsim_options);

    if p.SimulateMPCS
        fprintf_internal('\nSimulating MPCs...\n')
        mpc_simulator.solve(stats.pmf);
    end

    for ii = 1:6
        stats.other.sim_mpcs(ii).avg_0_quarterly = mpc_simulator.sim_mpcs(ii).quarterly;
        stats.other.sim_mpcs(ii).avg_0_annual = mpc_simulator.sim_mpcs(ii).annual;
    end

    clear mpc_simulator

    %% ----------------------------------------------------------------
    % SIMULATE MPCs OUT OF NEWS
    % -----------------------------------------------------------------
    shocks = [4,5,6];
    shockperiod = 4;
    mpc_simulator = HACTLib.computation.MPCSimulator(...
        p, income, grdKFE, KFE, shocks, shockperiod,...
        trans_dyn_solver.savedTimesUntilShock, p.mpcsim_options);

    if (p.SimulateMPCS_news == 1)
        fprintf_internal('\nSimulating MPCs...\n')
        mpc_simulator.solve(stats.pmf);
    end

    for ii = 1:6
        stats.other.sim_mpcs(ii).avg_4_quarterly = mpc_simulator.sim_mpcs(ii).quarterly;
        stats.other.sim_mpcs(ii).avg_4_annual = mpc_simulator.sim_mpcs(ii).annual;
    end

    clear mpc_simulator
    
    %% ----------------------------------------------------------------
    % MPCs WITHOUT INCOME RISK
    % -----------------------------------------------------------------
    mpc_finder_norisk = HACTLib.computation.MPCs(...
    	p, income_norisk, grdKFE_norisk, 'no_inc_risk', true);
    
    if (p.ComputeMPCS == 1) && (p.SolveNoRisk == 1)
    	fprintf_internal('\nComputing MPCs for model without income risk...\n')
    	mpc_finder_norisk.solve(KFE_nr, Au_nr);
    end
    stats.other.mpcs_nr = mpc_finder_norisk.mpcs;
    
    %% ----------------------------------------------------------------
    % DECOMPOSITIONS
    % -----------------------------------------------------------------
    import HACTLib.computation.Decomp
    decomp_obj = Decomp(p, grdKFE.b.vec, stats, income);

    if p.ComputeMPCS && p.OneAsset
        fprintf_internal('\nPerforming decompositions...\n')
        decomp_obj.compute();
    end

    stats.add_decomps(decomp_obj);
    
    %% ----------------------------------------------------------------
    % HOUSEKEEPING
    % -----------------------------------------------------------------
    import HACTLib.aux.to_structure

    grd.clean();
    grdKFE.clean();
    income.clean();

    grd = to_structure(grd);
    grdKFE = to_structure(grdKFE);
    p = to_structure(p);
    income = to_structure(income);

    if p.saveGrids
        save_grids(p, grd);
    end

	stats.clean();
    stats = HACTLib.aux.to_structure(stats);

    if final
        fname = sprintf('output_%d.mat', p.param_index);
        fpath = fullfile('output', fname);
        save(fpath,'stats','grd','grdKFE','p','KFE','income')
    end

    if save_iteration_results
        fname = sprintf('rho_ra_iterations_%d.mat', p.param_index);
        fpath = fullfile('output', fname);

        if p.calibrator.iter == 1
            results = struct();
            results.kappa1 = p.kappa1;
            results.kappa2 = p.kappa2;
            results.median_w = stats.median_totw.value;
            results.median_b = stats.median_liqw.value;
            results.rho = p.rho;
            results.r_a = p.r_a;
            results.rho_bounds = p.calibration_bounds{1};
            results.r_a_bounds = p.calibration_bounds{2};
        else
            results = load(fpath);
            results.kappa1(end+1) = p.kappa1;
            results.kappa2(end+1) = p.kappa2;
            results.median_w(end+1) = stats.median_totw.value;
            results.median_b(end+1) = stats.median_liqw.value;
            results.rho(end+1) = p.rho;
            results.r_a(end+1) = p.r_a;
        end

        save(fpath, '-struct', 'results')
    end

    clear solver.two_asset.solver
end

function save_grids(p, grd)
    nx = max(p.nb, p.na);
    bgrid = [grd.b.vec; NaN(nx-p.nb, 1)];
    agrid = [grd.a.vec; NaN(nx-p.na, 1)];
    grids = table(bgrid, agrid);
    fpath = fullfile('output',...
        sprintf('grids%d.xlsx', p.param_index));
    writetable(grids, fpath, 'WriteVariableNames', true);
end