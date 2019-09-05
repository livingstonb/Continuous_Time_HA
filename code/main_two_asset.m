function stats = main_two_asset(runopts)
    % Main function file for this repository. If IterateRho = 1, this script
    % first tries to find valid lower and upper bounds for rho (this may fail
    % in some cases), and then iterates on rho once valid bounds are found.

    %% --------------------------------------------------------------------
    % GET PARAMETERS AND CREATE GRID, INCOME OBJECTS
    % ---------------------------------------------------------------------
    if strcmp(runopts.mode,'grid_test')
    	p = setup.two_asset.params.grid_test_params(runopts);
    elseif strcmp(runopts.mode,'chi0_tests')
        p = setup.two_asset.params.chi0_tests(runopts);
    elseif strcmp(runopts.mode,'chi1_chi2_tests')
        p = setup.two_asset.params.chi1_chi2_tests(runopts);
    elseif strcmp(runopts.mode,'table_tests')
        p = setup.two_asset.params.table_tests(runopts);
    elseif strcmp(runopts.mode,'table_tests_bequests')
        p = setup.two_asset.params.table_tests_bequests(runopts);
    elseif strcmp(runopts.mode,'get_params')
		p = setup.two_asset.params.get_params(runopts);
    end
    p.print();
	
	dimsHJB = [p.nb p.na p.nz];
	dimsKFE = [p.nb_KFE p.na_KFE p.nz];
	income = setup.Income(runopts,p,dimsHJB,dimsKFE,false);
    income_norisk = setup.Income(runopts,p,dimsHJB,dimsKFE,true);

    p.update_ny(income.ny);

	grd = setup.two_asset.GridTwoAsset(p,income.ny,'HJB'); % grid for HJB
    grd_norisk = setup.two_asset.GridTwoAsset(p,1,'HJB');
	grdKFE = setup.two_asset.GridTwoAsset(p,income.ny,'KFE');% grid for KFE
    grdKFE_norisk = setup.two_asset.GridTwoAsset(p,1,'KFE');
    
    % check that borrowing limit does not violate NBL
    NBL = - min((1-p.wagetax-p.directdeposit)*income.y.vec+p.transfer) ...
        / (p.r_b_borr + p.deathrate*p.perfectannuities);
    msg = sprintf('bmin < natural borrowing limit (%f)',NBL);
    assert(p.bmin > NBL,msg);
    
	if p.IterateRho == 1
        % Look for valid rho lower bound and rho upper bound prior to AY iteration

        model_solver = @(x,y) solver.two_asset.solver(x,y,income,grd,grdKFE);
        aux.searchForRhoBounds
        
        %% ----------------------------------------------------------------
        % ITERATE OVER RHO TO MATCH MEAN ASSETS
        % -----------------------------------------------------------------
        runopts.RunMode = 'Iterate';
        iterate_rho = @(x) solver.two_asset.solver(runopts,p.reset_rho(x),income,grd,grdKFE);

        check_evals = @(x,y,z) aux.fzero_check(x,y,z,p);
        options = optimset('TolX',p.crit_AY,'OutputFcn',check_evals);
        [rho_final,~,exitflag] = fzero(iterate_rho,[rho_lb,rho_ub],options);

        if exitflag ~= 1
            error(['fzero failed, exitflag = ',num2str(exitflag)])
        end

        fprintf('\nIteration over rho completed.\n\n')
	else
        % Don't iterate
        rho_final = p.rho;
    end
    
    % Solve model one more time to get other output variables
    p.reset_rho(rho_final);
    runopts.RunMode = 'Final';
	[~,HJB,KFE,Au] = solver.two_asset.solver(runopts,p,income,grd,grdKFE);
    runopts.RunMode = 'NoRisk';
    [~,HJB_nr,KFE_nr,Au_nr] = solver.two_asset.solver(runopts,p,income_norisk,...
                                        grd_norisk,grdKFE_norisk);

    %% ----------------------------------------------------------------
    % COMPUTE STATISTICS
    % -----------------------------------------------------------------
    fprintf('\nComputing statistics\n')
    stats = statistics.two_asset.compute_statistics(p,income,grd,grdKFE,KFE);
    

    %% ----------------------------------------------------------------
    % COMPUTE MPCs
    % -----------------------------------------------------------------
    stats.mpcs = struct();
    dim2Identity = 'a';
    mpc_finder = statistics.MPCFinder(p,income,grdKFE,dim2Identity);
    trans_dyn_solver = solver.two_asset.TransitionalDynSolverTwoAsset(p,income,grdKFE);
    
    if p.ComputeMPCS == 1
    	fprintf('\nComputing MPCs out of an immediate shock...\n')
    	mpc_finder.solve(KFE,stats.pmf,Au);
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
    nperiods = 1;
    dim2Identity = 'a';
    mpc_simulator = statistics.two_asset.MPCSimulatorTwoAsset(...
    	p,income,grdKFE,KFE,shocks,nperiods,dim2Identity);

    if p.SimulateMPCS == 1
        fprintf('\nSimulating MPCs...\n')
        mpc_simulator.solve(p,income,grdKFE,stats.pmf);
    end
    for ii = 1:6
        stats.sim_mpcs(ii).avg_0_quarterly = mpc_simulator_immediateshock.sim_mpcs(ii).avg_quarterly;
        stats.sim_mpcs(ii).avg_0_annual = mpc_simulator_immediateshock.sim_mpcs(ii).avg_annual;
    end

    clear mpc_simulator
    
    %% ----------------------------------------------------------------
    % MPCs WITHOUT INCOME RISK
    % -----------------------------------------------------------------
    dim2Identity = 'a';
    mpc_finder_norisk = statistics.MPCFinder(...
    	p,income_norisk,grdKFE_norisk,dim2Identity);
    
    if p.ComputeMPCS == 1
    	fprintf('\nComputing MPCs for model without income risk...\n')
    	mpc_finder_norisk.solve(KFE_nr,stats.pmf_norisk,Au_nr);
    end
    stats.mpcs_nr = mpc_finder_norisk.mpcs;
    
    %% ----------------------------------------------------------------
    % DECOMPOSITIONS
    % -----------------------------------------------------------------
    fprintf('\nPerforming decompositions (if applicable)...\n')
    stats.decomp_norisk = statistics.two_asset.decomp_wrt_norisk(p,grdKFE,stats,income);
    stats.decompRA = statistics.two_asset.decompRA(p,grdKFE,stats);
    
    %% ----------------------------------------------------------------
    % HOUSEKEEPING
    % -----------------------------------------------------------------
    grd = aux.to_structure(grd);
    grdKFE = aux.to_structure(grdKFE);
    p = aux.to_structure(p);
    income = aux.to_structure(income);
    
    % empty some variables to reduce file size
    if ~strcmp(p.name,'baseline')
        grd.a.matrix = [];
        grd.b.matrix = [];
    end
    
    if runopts.ComputeMPCS == 1
	    stats.mpcs(1).mpcs_0_t = [];
	    stats.mpcs(2).mpcs_0_t = [];
	    stats.mpcs(3).mpcs_0_t = [];
	    stats.mpcs(4).mpcs_0_t = [];
	    stats.mpcs(6).mpcs_0_t = [];
	end
    
    save([runopts.savedir 'output_' runopts.suffix '.mat'],'stats','grd','grdKFE','p','KFE','income')
    clear solver.two_asset.solver
    fprintf('\nCode finished for the parameterization: \n\t%s\n\n',p.name)
end
