function decomp = decomp_wrt_norisk(p,grdKFE,stats,income)
    % Decomposition w.r.t. the no-risk model

    assets = grdKFE.a.matrix + grdKFE.b.matrix;
    assets = reshape(assets,p.nb_KFE*p.na_KFE,income.ny*p.nz);
    assets = assets(:,1);

    % initialize to NaN
    for ia = 1:numel(p.decomp_thresholds)
        decomp.term1(ia) = NaN;
        decomp.term2(ia) = NaN;
        decomp.term3(ia) = NaN;
        decomp.term4(ia) = NaN;
    end

    % check if required MPCs are available
    if (p.ComputeMPCS == 0) || (p.OneAsset == 0)
        return
    end
    
    m_ra = (p.rho + p.deathrate - p.deathrate*p.perfectannuities - p.r_b) / p.riskaver + p.r_b;
    m1 = stats.mpcs(5).mpcs(:,1);
    g1 = reshape(stats.pmf,[],income.ny*p.nz);

    % P(a & b)
    Pab = sum(reshape(g1,[],income.ny*p.nz),2);

    m1_wide = reshape(m1,[],income.ny*p.nz);
    m1_ab = m1_wide .* g1 ./ Pab;
    m1_ab = sum(m1_ab,2);
    % use arithmetic mean when denom < 1e-9
    Pab_small = Pab < 1e-9;
    m1_ab(Pab_small) = mean(m1_wide(Pab_small,:),2);

    % norisk model
    mbc = stats.mpcs_nr(5).mpcs(:,1);
    if sum(isnan(mbc(:))) + sum(isnan(m1(:))) > 0
        % contains NaNs
        return
    end

    % sum over income
    g1 = sum(g1,2);
    
    for ia = 1:numel(p.decomp_thresholds)
        zidx = assets <= p.decomp_thresholds(ia);
        
        decomp.term1(ia) = m_ra;
        decomp.term2(ia) = (m1_ab(zidx) - m_ra)' * g1(zidx);
        decomp.term3(ia) = (mbc(~zidx) - m_ra)' * g1(~zidx);
        decomp.term4(ia) = (m1_ab(~zidx) - mbc(~zidx))' * g1(~zidx);
    
    end
end