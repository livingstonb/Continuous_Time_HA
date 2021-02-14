function [rho_bds, r_a_bds] = get_rho_ra_bounds(kappa1, kappa2, printbds)
    if nargin < 3
        printbds = false;
    end

    r_b = 0.005;

    if (kappa2 == 0.1)
        if kappa1 <= 0.2
            rho_bds = [0.002, 0.007];
            r_a_bds = [0.006, 0.01];
        elseif kappa1 <= 0.5
            rho_bds = [0.005, 0.01];
            r_a_bds = [0.008, 0.015];
        elseif kappa1 <= 1.5
            rho_bds = [0.007, 0.03];
            r_a_bds = [0.01, 0.03];
        else
            rho_bds = [0.01, 0.03];
            r_a_bds = [0.03, 0.06];
        end
    elseif (kappa2 == 0.5)
        if kappa1 <= 0.2
            rho_bds = [0.002, 0.007];
            r_a_bds = [0.0055, 0.009];
        elseif kappa1 <= 0.5
            rho_bds = [0.003, 0.007];
            r_a_bds = [0.007, 0.015];
        elseif kappa1 <= 1.5
            rho_bds = [0.005, 0.01];
            r_a_bds = [0.0085, 0.015];
        else
            rho_bds = [0.0065, 0.02];
            r_a_bds = [0.011, 0.04];
        end
    elseif (kappa2 == 0.75)
        if kappa1 <= 0.2
            rho_bds = [0.0019, 0.006];
            r_a_bds = [0.0053, 0.009];
        elseif kappa1 <= 0.5
            rho_bds = [0.0021, 0.009];
            r_a_bds = [0.0053, 0.01];
        elseif kappa1 <= 1.5
            rho_bds = [0.0023, 0.01];
            r_a_bds = [0.006, 0.03];
        else
            rho_bds = [0.006, 0.02];
            r_a_bds = [0.011, 0.03];
        end
    elseif (kappa2 == 1)
        if kappa1 <= 1
            rho_bds = [0.0018, 0.004];
            r_a_bds = [0.0053, 0.0085];
        elseif kappa1 <= 2
            rho_bds = [0.003, 0.0065];
            r_a_bds = [0.0065, 0.01];
        elseif kappa1 <= 5
            rho_bds = [0.0045, 0.008];
            r_a_bds = [0.009, 0.012];
        else
            rho_bds = [0.006, 0.011];
            r_a_bds = [0.01, 0.014];
        end
    elseif (kappa2 == 2.0)
        if kappa1 <= 1
            rho_bds = [0.0018, 0.003];
            r_a_bds = [0.0052, 0.006];
        elseif kappa1 <= 2
            rho_bds = [0.0019, 0.003];
            r_a_bds = [0.0052, 0.0061];
        elseif kappa1 <= 5
            rho_bds = [0.0025, 0.004];
            r_a_bds = [0.0055, 0.007];
        else
            rho_bds = [0.0025, 0.004];
            r_a_bds = [0.0055, 0.007];
        end
    end

    if printbds
        fprintf('kappa1 = %f\n', kappa1)
        fprintf('kappa2 = %f\n', kappa2)
        fprintf('\trho_bds = [%f, %f]\n', rho_bds(1), rho_bds(2))
        fprintf('\tr_a_bds = [%f, %f]\n', r_a_bds(1), r_a_bds(2))
    end

    if rho_bds(1) > rho_bds(2)
        error('rho lb > rho ub')
    elseif rho_bds(1) < 0.001
        error('rho lb too small')
    elseif rho_bds(2) > 0.06
        error('rho ub too large')
    end

    if r_a_bds(1) > r_a_bds(2)
        error('r_a lb > r_a ub')
    elseif r_a_bds(1) < 0.00501
        error('r_a lb too small')
    elseif r_a_bds(2) > 0.06
        error('r_a ub too large')
    end
end