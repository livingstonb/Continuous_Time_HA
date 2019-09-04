%% ---------------------------------------------------------------
% SEARCH FOR VALID RHO LOWER BOUND
% ----------------------------------------------------------------
% Look for rho lower bound s.t. model converges and satisfies 
% assets > 3.5 in stationary distribution
disp('Looking for appropriate rho lower bound')

runopts.RunMode = 'Find rho_lb'; % Save V and g for first iterations
neg_stepsize = -1e-3;
pos_stepsize = 1e-3;
rho0 = p.rhoL;
rho_lb_finder = setup.RhoBoundsFinder(neg_stepsize,pos_stepsize,rho0);
n_bound_max = 40;
rhoLowerBoundFound = false;

p.reset_rho(rho0);

ii = 1;
while (ii <= n_bound_max) && (~rhoLowerBoundFound)
    try
        % Attempt to solve model
        AYdiff = model_solver(runopts,p,income,grd,grdKFE); 

        % If A/Y - target > 0, rho lb has been found
        if AYdiff > 0
            rhoLowerBoundFound = true;
        else
            rho_lb_finder.decrease(); % A/Y too low
        end              
    catch ME
        getReport(ME)
        rho_lb_finder.increase(); % A/Y probably too high
    end

    % update rho in parameters
    rho_lb = rho_lb_finder.get_new_rho();
    p.reset_rho(rho_lb);

    ii = ii + 1;
end

if ~rhoLowerBoundFound
    error('Rho lower bound not found')
end

%% ---------------------------------------------------------------
% SEARCH FOR VALID RHO UPPER BOUND
% ----------------------------------------------------------------
% Look for rho upper bound s.t. model converges and satisfies 
% assets < 3.5 in stationary distribution
disp('Looking for appropriate rho upper bound')

runopts.RunMode = 'Find rho_ub'; % Save V and g for first iterations
neg_stepsize = -2e-3;
pos_stepsize = 5e-3;
rhoH = rho_lb + 1e-3;
rho_ub_finder = setup.RhoBoundsFinder(neg_stepsize,pos_stepsize,rhoH);
n_bound_max = 40;
rhoUpperBoundFound = false;

p.reset_rho(rhoH);

ii = 1;
while (ii <= n_bound_max) && (~rhoUpperBoundFound)
    try
        % Attempt to solve model
        AYdiff = model_solver(runopts,p,income,grd,grdKFE); 

        % If A/Y - target < 0, rho ub has been found
        if AYdiff < 0
            rhoUpperBoundFound = true;
        else
            rho_ub_finder.increase(); % A/Y too high
        end              
    catch ME
        getReport(ME)
        rho_ub_finder.decrease(); % A/Y probably too low
    end

    % update rho in parameters
    rho_ub = rho_ub_finder.get_new_rho();
    p.reset_rho(rho_ub);

    ii = ii + 1;
end

if ~rhoUpperBoundFound
    error('Rho upper bound not found')
end

%% ----------------------------------------------------------------
% VALID RHO BOUNDS FOUND, NOW ITERATE OVER RHO TO MATCH MEAN ASSETS
% -----------------------------------------------------------------

runopts.RunMode = 'Iterate';
iterate_rho = @(x) solver.solver(runopts,p.reset_rho(x),income,grd,grdKFE);

check_evals = @(x,y,z) aux.fzero_check(x,y,z,p);
options = optimset('TolX',p.crit_AY,'OutputFcn',check_evals);
[rho_final,~,exitflag] = fzero(iterate_rho,[rho_lb,rho_ub],options);

if exitflag ~= 1
    error(['fzero failed, exitflag = ',num2str(exitflag)])
end

fprintf('\nIteration over rho completed.\n\n')