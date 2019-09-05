function decomp = decompRA(p,grdKFE,stats)
    % Decomposition w.r.t. the representative agent model

    assets = grdKFE.a.matrix(:) + grdKFE.b.matrix(:);

    % initialize to NaN
    decomp.RAmpc = NaN;
    decomp.Em1_less_mRA = NaN;
    decomp.term1 = NaN;
    decomp.term2 = NaN;
    decomp.term3 = NaN;

    % check if required MPCs are available
    if (p.ComputeMPCS == 0) || (p.OneAsset == 0)
        return
    end

    m1 = stats.mpcs(5).mpcs(:,1);
    pmf = stats.pmf(:);

    if sum(isnan(m1)) > 0
        return
    end

    % Decomp with respect to RA model
    m0 = (p.rho + p.deathrate - p.deathrate*p.perfectannuities - p.r_b) / p.riskaver + p.r_b;
    
    %% Find E[mpc|a=3.5]
    a_unique = unique(assets);
    Empc_a = zeros(size(a_unique));
    Psmall = false(size(a_unique));
    for ia = 1:numel(a_unique)
        a = a_unique(ia);
        Pa = sum(pmf(assets==a));
        % E[MPC|a]
        if Pa < 1e-9
            Psmall(ia) = true;
            Empc_a(ia) = NaN;
        else
            Empc_a(ia) = m1(assets==a)' * pmf(assets==a) / Pa;
        end
    end
    
    % ignore Empc_a for small values of Pa
    Empc_a = Empc_a(~Psmall);
    a_unique = a_unique(~Psmall);

    wsort = sortrows([a_unique Empc_a]);
    assets_sort = wsort(:,1);
    Empc_sort = wsort(:,2);

    winterp = griddedInterpolant(assets_sort,Empc_sort,'linear');
    mpc_atmean = winterp(3.5);

    decomp.RAmpc = m0;
    decomp.Em1_less_mRA = pmf' * m1 - m0;
    decomp.term1 = mpc_atmean - m0;
    decomp.term2 = 0;
    decomp.term3 = (m1 - m0)' * pmf - (mpc_atmean - m0);

end