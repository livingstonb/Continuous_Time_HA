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

    % get rid of z-dimension
    Em1 = stats1.mpcs(5).avg_0_quarterly(1);
    m1 = reshape(stats1.mpcs(5).mpcs(:,1),[p1.nb_KFE p1.na_KFE p1.nz income0.ny]);
    pmf1_z = stats1.pmf;
    % find P(z|b,a,y)
    P_bay = reshape(sum(pmf1_z,3),[p1.nb_KFE p1.na_KFE 1 income0.ny]);
    Pz_bay = pmf1_z ./ P_bay;
    m1 = sum(Pz_bay .* m1,3); % integrate out z
    m1(squeeze(P_bay)<1e-8) = 0; % avoid NaN
    m1 = m1(:);
    pmf1 = sum(pmf1_z,3);
    pmf1 = pmf1(:);

    Em0 = stats0.mpcs(5).avg_0_quarterly(1);
    m0 = stats0.mpcs(5).mpcs(:,1); % assume nz = 1 for baseline
    pmf0 = stats0.pmf(:);

    if all(isnan(m1)) || all(isnan(m0))
        return
    end

    % Main decomposition
    decomp.Em1_less_Em0 = Em1 - Em0;
    decomp.term1 = (m1 - m0)' * pmf0; 
    decomp.term2 = m0' * (pmf1 - pmf0);
    decomp.term3 = (m1 - m0)' * (pmf1 - pmf0);

    % get interpolant from assets to int_0^{epsilon}
    if p0.OneAsset == 1
	    m0g0interp = aux.interpolate_integral(assets0, m0, pmf0);
	    m1g0interp = aux.interpolate_integral(assets0, m1, pmf0);
	    m0g1interp = aux.interpolate_integral(assets0, m0, pmf1);
	    m1g1interp = aux.interpolate_integral(assets0, m1, pmf1);
	end

    % Decomposition of distribution effect
    for ia = 1:numel(p0.decomp_thresholds)
    	abar = p0.decomp_thresholds(ia); % HtM threshold

        if p0.OneAsset == 1
        	if abar == 0
	        	idx = assets0 <= abar;
	        	% HtM households
	        	decomp.term2a(ia) = m0(idx)' * (pmf1(idx) - pmf0(idx));
	        	% Non-HtM households
	        	decomp.term2b(ia) = m0(~idx)' * (pmf1(~idx) - pmf0(~idx));
	            
	            decomp.term2c(ia) = NaN;
	        else
	        	decomp.term2a(ia) = m0g1interp(abar) - m0g0interp(abar);
	        	decomp.term2b(ia) = (m0(:)' * pmf1(:) - m0g1interp(abar)) ...
	        		- (Em0 - m0g0interp(abar));
        	end
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