function [stats,p] = main(runopts, p)
    % Instantiates necessary classes and calls functions to solve the
    % model and compute statistics
    %
    % Parameters
    % ----------
    % runopts : a structure containing run options
    %
    % p : a Params object containing model parameters
    %
    % Returns
    % -------
    % stats : a structure containing statistics from the solved model
    %
    % p : the Params object used to solve the model

    import HACTLib.model_objects.Grid
    import HACTLib.model_objects.Income

    %% --------------------------------------------------------------------
    % CREATE GRID, INCOME OBJECTS
    % ---------------------------------------------------------------------
	income_path = fullfile(runopts.direc, 'input', p.income_dir);

    % Main income process
	income = Income(income_path, p, false);

    % Turn off income risk (set y equal to the mean)
    income_norisk = Income(runopts.direc, p, true);

    p.set("ny", income.ny, true);

    % Natural borrowing limit
    NBL = - min((1-p.wagetax-p.directdeposit)*income.y.vec+p.transfer) ...
            / (p.r_b_borr + p.deathrate*p.perfectannuities);
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

    if numel(p.rhos) > 1
    	grd.add_zgrid(p.rhos', p.na);
    	grd_norisk.add_zgrid(p.rhos', p.na);
    	grdKFE.add_zgrid(p.rhos', p.na_KFE);
    	grdKFE_norisk.add_zgrid(p.rhos', p.na_KFE);
    end

    % Add net income variables
    income.set_net_income(p, grd, grdKFE);
    income_norisk.set_net_income(p, grd_norisk, grdKFE_norisk);

    runopts.RunMode = 'Final';
    model = HACTLib.model_objects.Model(p, grd, grdKFE, income);
    model.initialize();
    [HJB, KFE, Au] = model.solve();

    if p.NoRisk == 1
        runopts.RunMode = 'NoRisk';
        model_nr = HACTLib.model_objects.Model(...
            p, grd_norisk, grdKFE_norisk, income_norisk);
        model_nr.initialize();
        [HJB_nr, KFE_nr, Au_nr] = model_nr.solve();
    end

    %% ----------------------------------------------------------------
    % COMPUTE STATISTICS
    % -----------------------------------------------------------------
    fprintf('\nComputing statistics\n')
    stats = statistics.compute_statistics(p,income,grd,grdKFE,KFE);

    %% ----------------------------------------------------------------
    % COMPUTE MPCs
    % -----------------------------------------------------------------
    stats.mpcs = struct();

    import HACTLib.computation.MPCs
    mpc_finder = MPCs(p, income, grdKFE, p.mpc_options);

    shocks = [4,5,6];
    import HACTLib.computation.MPCsNews
    trans_dyn_solver = MPCsNews(p, income, grdKFE, shocks, p.mpcs_news_options);
    
    if p.ComputeMPCS == 1
    	fprintf('\nComputing MPCs out of an immediate shock...\n')
        mpc_finder.solve(KFE, stats.pmf, Au);
    end
    for ishock = 1:6
        stats.mpcs(ishock).avg_0_quarterly = mpc_finder.mpcs(ishock).quarterly;
        stats.mpcs(ishock).avg_0_annual = mpc_finder.mpcs(ishock).annual;
        stats.mpcs(ishock).mpcs = mpc_finder.mpcs(ishock).mpcs;
    end
    
    if p.ComputeMPCS_news == 1
    	fprintf('Computing MPCs out of news...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf,mpc_finder.cum_con_baseline);
    elseif p.SimulateMPCS_news == 1
        fprintf('Iterating backward to find policy functions for future shock...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf);
    end
    for ishock = 1:6
        stats.mpcs(ishock).avg_1_quarterly = trans_dyn_solver.mpcs(ishock).avg_1_quarterly;
    	stats.mpcs(ishock).avg_4_quarterly = trans_dyn_solver.mpcs(ishock).avg_4_quarterly;
        stats.mpcs(ishock).avg_4_annual = trans_dyn_solver.mpcs(ishock).avg_4_annual;
    end


    %% ----------------------------------------------------------------
    % SIMULATE MPCs
    % -----------------------------------------------------------------
    stats.sim_mpcs = struct();
    shocks = [4,5,6];
    shockperiod = 0;
    mpc_simulator = HACTLib.computation.MPCSimulator(...
    	p, income, grdKFE, KFE, shocks, shockperiod, p.mpcsim_options);

    if p.SimulateMPCS == 1
        fprintf('\nSimulating MPCs...\n')
        mpc_simulator.solve(stats.pmf);
    end
    for ii = 1:6
        stats.sim_mpcs(ii).avg_0_quarterly = mpc_simulator.sim_mpcs(ii).avg_quarterly;
        stats.sim_mpcs(ii).avg_0_annual = mpc_simulator.sim_mpcs(ii).avg_annual;
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
        stats.sim_mpcs(ii).avg_4_quarterly = mpc_simulator.sim_mpcs(ii).avg_quarterly;
        stats.sim_mpcs(ii).avg_4_annual = mpc_simulator.sim_mpcs(ii).avg_annual;
    end

    clear mpc_simulator
    
    %% ----------------------------------------------------------------
    % MPCs WITHOUT INCOME RISK
    % -----------------------------------------------------------------
    mpc_finder_norisk = HACTLib.computation.MPCs(...
    	p, income_norisk, grdKFE_norisk);
    
    if p.ComputeMPCS == 1
    	fprintf('\nComputing MPCs for model without income risk...\n')
    	mpc_finder_norisk.solve(KFE_nr, stats.pmf_norisk, Au_nr);
    end
    stats.mpcs_nr = mpc_finder_norisk.mpcs;
    
    %% ----------------------------------------------------------------
    % DECOMPOSITIONS
    % -----------------------------------------------------------------
    fprintf('\nPerforming decompositions (if applicable)...\n')
    stats.decomp_norisk = statistics.decomp_wrt_norisk(p, grdKFE, stats, income);
    stats.decompRA = statistics.decompRA(p, grdKFE, stats);
    
    %% ----------------------------------------------------------------
    % HOUSEKEEPING
    % -----------------------------------------------------------------
    import HACTLib.aux.to_structure
    grd = to_structure(grd);
    grdKFE = to_structure(grdKFE);
    p = to_structure(p);
    income = to_structure(income);
    
    % empty some variables to reduce file size
    if ~strcmp(p.name,'baseline')
        grd.a.matrix = [];
        grd.b.matrix = [];
    end
    
    if p.ComputeMPCS
	    stats.mpcs(1).mpcs_0_t = [];
	    stats.mpcs(2).mpcs_0_t = [];
	    stats.mpcs(3).mpcs_0_t = [];
	    stats.mpcs(4).mpcs_0_t = [];
	    stats.mpcs(6).mpcs_0_t = [];
	end
    
    if p.SaveResults == 1
    	spath = fullfile(runopts.savedir, ['output_' runopts.suffix '.mat']);
        save(spath,'stats','grd','grdKFE','p','KFE','income')
    end
    clear solver.two_asset.solver
    fprintf('\nCode finished for the parameterization: \n\t%s\n\n',p.name)
end
