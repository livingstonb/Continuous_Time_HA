function [stats, stats_alt] = main(p, save_results)
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

    if p.OneAsset
        p.set("ComputeMPCS_illiquid", false, true);
    end

    %% --------------------------------------------------------------------
    % CREATE GRID, INCOME OBJECTS
    % ---------------------------------------------------------------------
	income_path = fullfile('input', p.income_dir);

    % Main income process
	income = Income(income_path, p, false);

    % Turn off income risk (set y equal to the mean)
    income_norisk = Income('', p, true);

    p.set("ny", income.ny, true);

    % Natural borrowing limit
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

    model = HACTLib.model_objects.Model(p, grd, grdKFE, income);
    model.initialize();
    [~, KFE, Au] = model.solve();

    if p.NoRisk == 1
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
    fprintf('\nComputing statistics\n')    
    stats = HACTLib.Statistics(p, income, grdKFE, KFE);
    
    kernel_options.ktype = 'gaussian';
    kernel_options.h = 0.2;
    kernel_options.force_fit_cdf_low = [];
    kernel_options.rescale_and_log = true;
    stats.compute_statistics(kernel_options);

    %% ----------------------------------------------------------------
    % COMPUTE MPCs
    % -----------------------------------------------------------------
    import HACTLib.computation.MPCs

    mpc_finder = MPCs(p, income, grdKFE, p.mpc_options);
    mpc_finder_illiquid = MPCs(p, income, grdKFE, p.mpc_options_illiquid);

    shocks = [4, 5, 6];
    import HACTLib.computation.MPCsNews
    trans_dyn_solver = MPCsNews(p, income, grdKFE, shocks, p.mpcs_news_options);
    
    if p.ComputeMPCS
    	fprintf('\nComputing MPCs out of an immediate shock...\n')
        mpc_finder.solve(KFE, Au, stats.pmf);
    end

    if p.ComputeMPCS_illiquid
        fprintf('\nComputing illiquid MPCs out of an immediate shock...\n')
        mpc_finder_illiquid.solve(KFE, Au, stats.pmf);
    end
    
    stats.add_mpcs(mpc_finder);
    stats.add_mpcs(mpc_finder_illiquid);
    
    if p.ComputeMPCS_news == 1
    	fprintf('Computing MPCs out of news...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf,mpc_finder.cum_con_baseline);
    elseif p.SimulateMPCS_news == 1
        fprintf('Iterating backward to find policy functions for future shock...\n')
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

    if p.SimulateMPCS == 1
        fprintf('\nSimulating MPCs...\n')
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

    if p.SimulateMPCS_news == 1
        fprintf('\nSimulating MPCs...\n')
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
    
    if (p.ComputeMPCS == 1) && (p.NoRisk == 1)
    	fprintf('\nComputing MPCs for model without income risk...\n')
    	mpc_finder_norisk.solve(KFE_nr, Au_nr);
    end
    stats.other.mpcs_nr = mpc_finder_norisk.mpcs;
    
    %% ----------------------------------------------------------------
    % DECOMPOSITIONS
    % -----------------------------------------------------------------
    import HACTLib.computation.Decomp
    decomp_obj = Decomp(p, grdKFE.b.vec, stats, income);

    if p.ComputeMPCS && p.OneAsset
        fprintf('\nPerforming decompositions...\n')
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
        nx = max(p.nb, p.na);
        bgrid = [grd.b.vec; NaN(nx-p.nb, 1)];
        agrid = [grd.a.vec; NaN(nx-p.na, 1)];
        grids = table(bgrid, agrid);
        fpath = fullfile('output',...
            sprintf('grids%d.xlsx', p.param_index));
        writetable(grids, fpath, 'WriteVariableNames', true);
    end

	stats.clean();
    stats = HACTLib.aux.to_structure(stats);

    if save_results
        fname = sprintf('output_%d.mat', p.param_index);
        fpath = fullfile('output', fname);
        save(fpath,'stats','grd','grdKFE','p','KFE','income')
    end

    clear solver.two_asset.solver
    fprintf('\nCode finished for the parameterization: \n\t%s\n\n',p.name)
end
