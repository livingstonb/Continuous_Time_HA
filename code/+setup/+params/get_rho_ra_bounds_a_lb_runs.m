function [rho_bds, r_a_bds] = get_rho_ra_bounds_a_lb_runs(a_lb, printbds)
    if nargin < 2
        printbds = false;
    end

    if a_lb <= 50000
        rho_bds = [0.06, 0.08];
        r_a_bds = [0.055, 0.08];
    else
        rho_bds = [];
        r_a_bds = [];
    end

    % [100, 500, 1000, 2000, 5000, 10000, 20000]

%     rho_bds = [0.008, 0.013];
%     r_a_bds = [0.06, 0.13];

    if printbds
        % fprintf('kappa1 = %f\n', kappa1)
        % fprintf('kappa2 = %f\n', kappa2)
        fprintf('a_lb = %f\n', a_lb)
        fprintf('\trho_bds = [%f, %f]\n', rho_bds(1), rho_bds(2))
        fprintf('\tr_a_bds = [%f, %f]\n', r_a_bds(1), r_a_bds(2))
    end

    % if rho_bds(1) > rho_bds(2)
    %     error('rho lb > rho ub')
    % elseif rho_bds(1) < 0.001
    %     error('rho lb too small')
    % elseif rho_bds(2) > 0.5
    %     error('rho ub too large')
    % end

    % if r_a_bds(1) > r_a_bds(2)
    %     error('r_a lb > r_a ub')
    % elseif r_a_bds(1) < 0.00501
    %     error('r_a lb too small')
    % elseif r_a_bds(2) > inf
    %     error('r_a ub too large')
    % end
%     rho_bds
%     r_a_bds
end