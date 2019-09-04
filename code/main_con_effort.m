function [stats,p,grdKFE,KFE] = main_con_effort(runopts)
    % Main function file for this repository. If IterateRho = 1, this script
    % first tries to find valid lower and upper bounds for rho (this may fail
    % in some cases), and then iterates on rho once valid bounds are found.

    %% --------------------------------------------------------------------
    % GET PARAMETERS AND CREATE GRID, INCOME OBJECTS
    % ---------------------------------------------------------------------

    p = setup.con_effort.get_params(runopts);
    p.print();
    
    if p.b_gcurv_neg < 1
    	warning('using a curved negative asset grid may violate min_grid_spacing')
    end
	
    dimsHJB = [p.nb p.nc p.ny];
    dimsKFE = [p.nb_KFE p.nc_KFE p.ny];
	income = setup.Income(runopts,p,dimsHJB,dimsKFE,false);
    p.update_ny(income.ny);
    
    dimsHJB = [p.nb p.nc 1];
    dimsKFE = [p.nb_KFE p.nc_KFE 1];
	income_norisk = setup.Income(runopts,p,dimsHJB,dimsKFE,true);

	grd = setup.con_effort.GridConEffort(p,p.ny,'HJB'); % grid for HJB
    grdKFE = setup.con_effort.GridConEffort(p,p.ny,'KFE'); % grid for KFE
    grd_norisk = setup.con_effort.GridConEffort(p,1,'HJB');
    grdKFE_norisk = setup.con_effort.GridConEffort(p,1,'KFE');
    
	if runopts.IterateRho == 1
        model_solver = @(x,y) solver.con_effort.solver(x,y,income,grd,grdKFE);
        aux.searchForRhoBounds.m
	else
        % Don't iterate
        rho_final = p.rho;
    end
    
    %% ----------------------------------------------------------------
    % FINAL RUN
    % -----------------------------------------------------------------
    p.reset_rho(rho_final);
    runopts.RunMode = 'Final';
    % update grd, grdKFE with c-grids
	[AYdiff,HJB,KFE,Au,grd,grdKFE] = solver.con_effort.solver(runopts,p,income,grd,grdKFE);
    
    %% --------------------------------------------------------------------
    % POLICY FUNCTION PLOTS
    % ---------------------------------------------------------------------
    % index of consumption at 0, 25, 50, 75, 100 pctildes of cgrid
    % (based on cgrid only, not equilibrium distribution of consumption)
    % qind = [2; round(nc/4); round(nc/2); round(nc*3/4); nc-1];
    
    % % h for lowest income
    % for iq = qind
    %     plot(grd.a.vec,HJB.h(:,iq,1))
    %     hold on
    % end
    % title('Lowest income')
    % legend('Cmin','(3/4)*Cmin + (1/4)*Cmax','(1/2)*Cmin + (1/2)*Cmax','(1/4)*Cmin + (3/4)*Cmax','Cmax')
    % xlabel('Assets, a')
    % ylabel(['h = ' p.hdef])
    % axis([0 10 -inf inf])
    % savepath = [runopts.savedir 'policy_ylow_param' runopts.suffix '.png'];
    % saveas(gcf,savepath);
    % close;
    
    % % h for highest income
    % for iq = qind
    %     plot(grd.a.vec,HJB.h(:,iq,27))
    %     hold on
    % end
    % title('Middle income')
    % legend('Cmin','(3/4)*Cmin + (1/4)*Cmax','(1/2)*Cmin + (1/2)*Cmax','(1/4)*Cmin + (3/4)*Cmax','Cmax')
    % xlabel('Assets, a')
    % ylabel(['h = ' p.hdef])
    % axis([0 10 -inf inf])
    % savepath = [runopts.savedir 'policy_ymed_param' runopts.suffix '.png'];
    % saveas(gcf,savepath);
    % close;
    
    % % h for highest income
    % for iq = qind
    %     plot(grd.a.vec,HJB.h(:,iq,55))
    %     hold on
    % end
    % title('Highest income')
    % legend('Cmin','(3/4)*Cmin + (1/4)*Cmax','(1/2)*Cmin + (1/2)*Cmax','(1/4)*Cmin + (3/4)*Cmax','Cmax')
    % xlabel('Assets, a')
    % ylabel(['h = ' p.hdef])
    % axis([0 10 -inf inf])
    % savepath = [runopts.savedir 'policy_yhigh_param' runopts.suffix '.png'];
    % saveas(gcf,savepath);
    % close;

    %% ----------------------------------------------------------------
    % COMPUTE STATISTICS
    % -----------------------------------------------------------------
    fprintf('Computing statistics\n\n')
    stats = statistics.con_effort.compute_statistics(p,income,grdKFE,KFE);

    %% ----------------------------------------------------------------
    % COMPUTE MPCs
    % -----------------------------------------------------------------
    mpc_finder = statistics.MPCFinder(p,income,grdKFE,'c');
    trans_dyn_solver = statistics.con_effort.TransitionalDynSolverConEffort(p,income,grdKFE);
    if p.ComputeMPCS == 1
    	fprintf('Computing MPCs out of an immediate shock...\n')
        mpc_finder.solve(KFE,stats.ptmass,Au);
    end

    if p.ComputeMPCS_news == 1
    	fprintf('Computing MPCs out of news...\n')
    	trans_dyn_solver.solve(KFE,stats.ptmass,mpc_finder.cum_con_baseline);
    elseif p.SimulateMPCS_news == 1
        fprintf('Iterating backward to find policy functions for future shock...\n')
    	trans_dyn_solver.solve(KFE,stats.ptmass);
    end

    stats.mpcs = mpc_finder.mpcs;
    for ii = 1:6
    	stats.mpcs(ii).avg_1_t = trans_dyn_solver.mpcs(ii).avg_1_t;
    	stats.mpcs(ii).avg_4_t = trans_dyn_solver.mpcs(ii).avg_4_t;
    end

    %% ----------------------------------------------------------------
    % SIMULATE MPCs
    % -----------------------------------------------------------------
    % mpcs out of immediate shock
    shocks = [4,5,6];
    nquarters = 4;
    mpc_simulator = statistics.MPCSimulatorImmediateShock(...
    	p,income,grdKFE,KFE,shocks,nquarters);
    if p.SimulateMPCS == 1
    	fprintf('\nSimulating MPCs...\n')
        mpc_simulator.solve(p,income,grdKFE,stats.ptmass);
    end

    stats.sim_mpcs = mpc_simulator.sim_mpcs;
    clear mpc_simulator;

    % mpcs out of news
    shocks = [5];
    
    nquarters = 4;
    shockperiod = 4;
    mpc_simulator_news_4 = statistics.MPCSimulatorFutureShock(...
    	p,income,grdKFE,KFE,shocks,nquarters,shockperiod);
    
    nquarters = 1;
    shockperiod = 1;
    mpc_simulator_news_1 = statistics.MPCSimulatorFutureShock(...
    	p,income,grdKFE,KFE,shocks,nquarters,shockperiod);
    
    if p.SimulateMPCS_news == 1
    	fprintf('\nSimulating MPCs out of news...\n')
        mpc_simulator_news_1.solve(p,income,grdKFE,stats.ptmass);
        mpc_simulator_news_4.solve(p,income,grdKFE,stats.ptmass);
    end
    for ii = 1:6
    	stats.sim_mpcs(ii).avg_1_t = mpc_simulator_news_1.sim_mpcs(ii).avg_1_t;
    	stats.sim_mpcs(ii).avg_4_t = mpc_simulator_news_4.sim_mpcs(ii).avg_4_t;
    end
    clear mpc_simulator_news_1 mpc_simulator_news_4

    %% ----------------------------------------------------------------
    % HOUSEKEEPING
    % -----------------------------------------------------------------
    grd = aux.to_structure(grd);
    grdKFE = aux.to_structure(grdKFE);
    p = aux.to_structure(p);
    
    % empty some variables to reduce file size
    if ~strcmp(p.name,'baseline')
        grd.a.matrix = [];
        grd.c.matrix = [];
    end
    
    grd.a.dF = [];
    grd.a.dB = [];
    grd.a.d_tilde = [];
    grdKFE.a.dF = [];
    grdKFE.a.dB = [];
    grdKFE.a.d_tilde = [];
    
    grd.c.dF = [];
    grd.c.dB = [];
    grd.c.d_tilde = [];
    grd.dac_tilde = [];
    grdKFE.c.dF = [];
    grdKFE.c.dB = [];
    grdKFE.c.d_tilde = [];
    grdKFE.dac_tilde = [];
    
    save([runopts.savedir 'output_' runopts.suffix '.mat'],...
        'stats','grd','grdKFE','p','KFE','income')
    clear solver.con_effort.solver
    
    %% ----------------------------------------------------------------
    % OUTPUT TO TERMINAL
    % -----------------------------------------------------------------
    if (p.ComputeMPCS==1) && (p.SimulateMPCS==1)
        fprintf('\nSome results...\n')
        for quarter = 1:4
            fprintf('    --- Quarter %i MPC out of 0.01 shock ---\n',quarter)
            fprintf('\tFeynman-Kac = %f\n',stats.mpcs(5).avg_0_t(quarter))
            fprintf('\tSimulated   = %f\n',stats.sim_mpcs(5).avg_0_t(quarter))
        end
    end
    
    fprintf('\nCode finished for the parameterization: \n\t%s\n\n',p.name)
end
