function decomp = decomp_baseline(s0,s1)
    % Decomposition of E[mpc1] - E[mpc0]
    % 0 is baseline
    % 1 is the experiment
    
    p0 = s0.p;
    stats0 = s0.stats;
    p1 = s1.p;
    stats1 = s1.stats;
    
    decomp.Em1_less_Em0 = NaN;
    decomp.term1 = NaN;
    decomp.term2 = NaN;
    decomp.term3 = NaN;
    decomp.term2a = NaN(numel(p0.decomp_thresholds),1);
    decomp.term2b = NaN(numel(p0.decomp_thresholds),1);
    decomp.term2c = NaN(numel(p0.decomp_thresholds),1);
    
    if isequaln(s0,s1)
        return
    end

    grdKFE0 = s0.grdKFE;
    
    income0 = s0.income;

    assets0 = grdKFE0.a.matrix(:) + grdKFE0.b.matrix(:);

    m1 = stats1.mpcs(5).mpcs(:,1);
    pmf1 = stats1.pmf(:);
    
    % get rid of z-dimension if present
    if p1.nz > 1
        dims = [p1.nb_KFE p1.na_KFE income0.ny p1.nz];
        m1 = reshape(m1,[],p1.nz);
        pmf1 = reshape(pmf1,[],p1.nz);
        
        m1 = sum(m1 .* pmf1,2) ./ sum(pmf1,2);
        m1 = m1(:);
        m1(sum(pmf1,2)<1e-7) = 0;
        pmf1 = sum(pmf1,2);
    end

    m0 = stats0.mpcs(5).mpcs(:,1);
    pmf0 = stats0.pmf(:);

    if all(isnan(m1)) || all(isnan(m0))
        return
    end

    % Main decomposition
    decomp.Em1_less_Em0 = stats1.mpcs(5).avg_0_t(1) - stats0.mpcs(5).avg_0_t(1);
    decomp.term1 = (m1 - m0)' * pmf0; 
    decomp.term2 = m0' * (pmf1 - pmf0);
    decomp.term3 = (m1 - m0)' * (pmf1 - pmf0);

    % Decomposition of distribution effect
    for ia = 1:numel(p0.decomp_thresholds)
    	abar = p0.decomp_thresholds(ia); % HtM threshold

        if p0.OneAsset == 1
        	idx = assets0 <= abar;
        	% HtM households
        	decomp.term2a(ia) = m0(idx)' * (pmf1(idx) - pmf0(idx));
        	% Non-HtM households
        	decomp.term2b(ia) = m0(~idx)' * (pmf1(~idx) - pmf0(~idx));
            
            decomp.term2c(ia) = NaN;
        else
            loc_zero = (grdKFE0.a.matrix(:)<=abar) & (grdKFE0.b.matrix(:)<=abar);
            decomp.term2a(ia) = m0(loc_zero)' * (pmf1(loc_zero) - pmf0(loc_zero));

            apos = grdKFE0.a.matrix(:) > abar;
            loc_bzero = grdKFE0.b.matrix(:) <= abar;
            idx = apos & loc_bzero;
            decomp.term2b(ia) = m0(idx)' * (pmf1(idx) - pmf0(idx));

            bpos = grdKFE0.b.matrix(:) > abar;
            decomp.term2c(ia) = m0(bpos)' * (pmf1(bpos) - pmf0(bpos));
        end
    end
end