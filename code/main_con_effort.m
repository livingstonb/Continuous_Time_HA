function [stats,p,grdKFE,KFE] = main_con_effort(runopts)
    % Main function file for this repository. If IterateRho = 1, this script
    % first tries to find valid lower and upper bounds for rho (this may fail
    % in some cases), and then iterates on rho once valid bounds are found.

    %% --------------------------------------------------------------------
    % GET PARAMETERS AND CREATE GRID, INCOME OBJECTS
    % ---------------------------------------------------------------------


    if strcmp(runopts.params_file,'get_params')
        p = setup.con_effort.get_params(runopts);
    elseif strcmp(runopts.params_file,'get_params2')
        p = setup.con_effort.get_params2(runopts);
    end 
    p.print();
    
    if p.b_gcurv_neg < 1
    	warning('using a curved negative asset grid may violate min_grid_spacing')
    end
	
    dimsHJB = [p.nb p.nc];
    dimsKFE = [p.nb_KFE p.nc_KFE];
	income = setup.Income(runopts,p,dimsHJB,dimsKFE,false);
    p.update_ny(income.ny);
    
    dimsHJB = [p.nb p.nc];
    dimsKFE = [p.nb_KFE p.nc_KFE];
	income_norisk = setup.Income(runopts,p,dimsHJB,dimsKFE,true);

	grd = setup.con_effort.GridConEffort(p,p.ny,'HJB'); % grid for HJB
    grdKFE = setup.con_effort.GridConEffort(p,p.ny,'KFE'); % grid for KFE
    grd_norisk = setup.con_effort.GridConEffort(p,1,'HJB');
    grdKFE_norisk = setup.con_effort.GridConEffort(p,1,'KFE');

    if numel(p.rhos) > 1
    	grd.add_zgrid(p.rhos',p.nc);
    	grd_norisk.add_zgrid(p.rhos',p.nc);
    	grdKFE.add_zgrid(p.rhos',p.nc_KFE);
    	grdKFE_norisk.add_zgrid(p.rhos',p.nc_KFE);
    end
    
	if runopts.IterateRho == 1
        model_solver = @(x,y) solver.con_effort.solver(x,y,income,grd,grdKFE);
        aux.searchForRhoBounds
        
        %% ----------------------------------------------------------------
        % VALID RHO BOUNDS FOUND, NOW ITERATE OVER RHO TO MATCH MEAN ASSETS
        % -----------------------------------------------------------------
        runopts.RunMode = 'Iterate';
        iterate_rho = @(x) solver.con_effort.solver(runopts,p.reset_rho(x),income,grd,grdKFE);

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
    
    %% ----------------------------------------------------------------
    % FINAL RUN
    % -----------------------------------------------------------------
    p.reset_rho(rho_final);
    runopts.RunMode = 'Final';
	[AYdiff,HJB,KFE,Au] = solver.con_effort.solver(runopts,p,income,grd,grdKFE);
    
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
    shocks = [2,3,5,6];
    trans_dyn_solver = solver.con_effort.TransitionalDynSolverConEffort(p,income,grdKFE,shocks);
    if p.ComputeMPCS == 1
    	fprintf('Computing MPCs out of an immediate shock...\n')
        mpc_finder.solve(KFE,stats.pmf,Au);
    end
    for ii = 1:6
        stats.mpcs(ii).avg_0_quarterly = mpc_finder.mpcs(ii).quarterly;
        stats.mpcs(ii).avg_0_annual = mpc_finder.mpcs(ii).annual;
        stats.mpcs(ii).mpcs = mpc_finder.mpcs(ii).mpcs;
    end

    if p.ComputeMPCS_news == 1
    	fprintf('Computing MPCs out of news...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf,mpc_finder.cum_con_baseline);
    elseif p.SimulateMPCS_news == 1
        fprintf('Iterating backward to find policy functions for future shock...\n')
    	trans_dyn_solver.solve(KFE,stats.pmf);
    end
    
    for ii = 1:6
    	stats.mpcs(ii).avg_1_quarterly = trans_dyn_solver.mpcs(ii).avg_1_quarterly;
    	stats.mpcs(ii).avg_4_quarterly = trans_dyn_solver.mpcs(ii).avg_4_quarterly;
        stats.mpcs(ii).avg_4_annual = trans_dyn_solver.mpcs(ii).avg_4_annual;
    end

    %% ----------------------------------------------------------------
    % SIMULATE MPCs
    % -----------------------------------------------------------------
    % mpcs out of immediate shock
    shocks = [2,3,5,6];
    mpc_simulator_immediateshock = statistics.con_effort.MPCSimulatorConEffort(...
    	p,income,grdKFE,KFE,shocks,0);
    if p.SimulateMPCS == 1
    	fprintf('\nSimulating MPCs...\n')
        mpc_simulator_immediateshock.solve(stats.pmf);
    end

    % mpcs out of news
    shocks = [5,6];
    mpc_simulator_q1shock = statistics.con_effort.MPCSimulatorConEffort(...
    	p,income,grdKFE,KFE,shocks,1);
    
    mpc_simulator_q4shock = statistics.con_effort.MPCSimulatorConEffort(...
    	p,income,grdKFE,KFE,shocks,4);
    
    if p.SimulateMPCS_news == 1
        mpc_simulator_q1shock.solve(stats.pmf);
        mpc_simulator_q4shock.solve(stats.pmf);
    end
    for ii = 1:6
    	stats.sim_mpcs(ii).responders_0_quarterly = mpc_simulator_immediateshock.sim_mpcs(ii).responders_quarterly;
        stats.sim_mpcs(ii).responders_0_annual= mpc_simulator_immediateshock.sim_mpcs(ii).responders_annual;
        stats.sim_mpcs(ii).avg_0_quarterly = mpc_simulator_immediateshock.sim_mpcs(ii).avg_quarterly;
        stats.sim_mpcs(ii).avg_0_quarterly_pos = mpc_simulator_immediateshock.sim_mpcs(ii).avg_quarterly_pos;
        stats.sim_mpcs(ii).avg_0_annual = mpc_simulator_immediateshock.sim_mpcs(ii).avg_annual;
        stats.sim_mpcs(ii).avg_0_annual_pos = mpc_simulator_immediateshock.sim_mpcs(ii).avg_annual_pos;
        stats.sim_mpcs(ii).responders_1_quarterly = mpc_simulator_q1shock.sim_mpcs(ii).responders_quarterly;
    	stats.sim_mpcs(ii).avg_1_quarterly = mpc_simulator_q1shock.sim_mpcs(ii).avg_quarterly;
    	stats.sim_mpcs(ii).avg_1_quarterly_pos = mpc_simulator_q1shock.sim_mpcs(ii).avg_quarterly_pos;
    	stats.sim_mpcs(ii).responders_4_quarterly = mpc_simulator_q4shock.sim_mpcs(ii).responders_quarterly;
        stats.sim_mpcs(ii).responders_4_annual = mpc_simulator_q4shock.sim_mpcs(ii).responders_annual;
    	stats.sim_mpcs(ii).avg_4_quarterly = mpc_simulator_q4shock.sim_mpcs(ii).avg_quarterly;
        stats.sim_mpcs(ii).avg_4_quarterly_pos = mpc_simulator_q4shock.sim_mpcs(ii).avg_quarterly_pos;
        stats.sim_mpcs(ii).avg_4_annual = mpc_simulator_q4shock.sim_mpcs(ii).avg_annual;
    end
    clear mpc_simulator_q1shock mpc_simulator_q4shock mpc_simulator_immediateshock

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
            fprintf('\tFeynman-Kac = %f\n',stats.mpcs(5).avg_0_quarterly(quarter))
            fprintf('\tSimulated   = %f\n',stats.sim_mpcs(5).avg_0_quarterly(quarter))
            fprintf('    --- Quarter %i MPC out of 0.1 shock ---\n',quarter)
            fprintf('\tFeynman-Kac = %f\n',stats.mpcs(6).avg_0_quarterly(quarter))
            fprintf('\tSimulated   = %f\n',stats.sim_mpcs(6).avg_0_quarterly(quarter))
        end
    end
    
    fprintf('\nCode finished for the parameterization: \n\t%s\n\n',p.name)
end
