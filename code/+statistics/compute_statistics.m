function stats = statistics(p, income, grd, grdKFE, KFE)
	% computes various statistics after solving the model

    import HACTLib.aux.compute_pct
    import HACTLib.aux.find_constrained
    import HACTLib.aux.direct_gini
	
    %% --------------------------------------------------------------------
    % Consumption by hand-to-mouth status
    % ---------------------------------------------------------------------
	na_KFE = p.na_KFE;
	nb_KFE = p.nb_KFE;
	ny = p.ny;
    nz = p.nz;

	stats.pmf = KFE.g .* grdKFE.trapezoidal.matrix;
    stats.pmf_norisk = sum(stats.pmf,4);
    stats.rho = p.rho;
    stats.beta_annualized = exp(-4*p.rho);
	wealth_mat = grdKFE.b.matrix + grdKFE.a.matrix;

    %% --------------------------------------------------------------------
    % WEALTH LEVELS
    % ---------------------------------------------------------------------
    stats.illiqw = stats.pmf(:)' * grdKFE.a.matrix(:);
    stats.liqw = stats.pmf(:)' * grdKFE.b.matrix(:);
    stats.totw = stats.illiqw + stats.liqw;
    
    %% --------------------------------------------------------------------
    % WEALTH PERCENTILES
    % ---------------------------------------------------------------------
    % g(b)
    tmp = reshape(stats.pmf, p.nb_KFE, []);
    stats.pmf_b = sum(tmp, 2);
    stats.cdf_b = cumsum(stats.pmf_b);
    
    % g(a)
    tmp = reshape(stats.pmf, p.nb_KFE, p.na_KFE, []);
    stats.pmf_a = squeeze(sum(sum(tmp, 3), 1));
    stats.cdf_a = cumsum(stats.pmf_a);
    
    % pmf(wealth)
    tmp = reshape(stats.pmf, p.nb_KFE*p.na_KFE, []);
    stats.pmf_wealth = reshape(sum(tmp, 2), [], p.na_KFE);
    wealth_values = reshape(wealth_mat(:,:,1,1), [], 1);

    tmp = sort([wealth_values, stats.pmf_wealth(:)]);
    [stats.values_cdf_wealth, iu] = unique(tmp(:,1), 'last');
    stats.cdf_wealth = cumsum(tmp(:,2));
    stats.cdf_wealth = stats.cdf_wealth(iu);
  
    % Interpolants for percentiles
    [cdf_b_u, iu] = unique(stats.cdf_b, 'first');
    lwinterp = griddedInterpolant(...
        cdf_b_u, grdKFE.b.vec(iu), 'pchip', 'nearest');

    if ~p.OneAsset
        [cdf_a_u, iu] = unique(stats.cdf_a, 'first');
        lwinterp = griddedInterpolant(...
            cdf_a_u, grdKFE.a.vec(iu), 'pchip', 'nearest');
    else
        iwinterp = @(a) NaN;
    end

    [cdf_w_u, iu] = unique(stats.cdf_wealth, 'first');
    winterp = griddedInterpolant(...
        cdf_w_u, stats.values_cdf_wealth(iu),..
        'pchip', 'nearest');

    pct_vals = p.wpercentiles / 100;
    stats.lwpercentile = lwinterp(pct_vals);
    stats.iwpercentile = iwinterp(pct_vals);
    stats.wpercentile = zinterp(pct_vals);

    stats.median_liqw = lwinterp(0.5);
    stats.median_illliqw = iwinterp(0.5);
    stats.median_totw = winterp(0.5);

    % Top shares
    wealth_values = reshape(wealth_mat(:,:,1,1), [], 1);
    tmp = sort([wealth_values, stats.pmf_wealth(:)]);
    % Amount of total assets that reside in each pt on sorted asset space
    totassets = tmp(:,1) .* tmp(:,2);
    % Fraction of total assets in each pt on asset space
    cumassets = cumsum(totassets) / stats.totw;
    [cumassets, iu] = unique(cumassets, 'last');
    cdf_u = cumsum(tmp(:,2));
    cdf_u = cdf_u(iu);

    cumwealthshare_interp = griddedInterpolant(...
        cdf_u, cumassets, 'pchip', 'nearest');
    stats.top10share = 1 - cumwealthshare_interp(0.9);
    stats.top1share = 1 - cumwealthshare_interp(0.99);

    %% --------------------------------------------------------------------
    % CONSTRAINED HOUSEHOLDS
    % ---------------------------------------------------------------------
    n_HtM = numel(p.epsilon_HtM);
    stats.constrained = zeros(1, n_HtM);
    stats.constrained_liq = zeros(1, n_HtM);

    % saving = 0
    temp = stats.pmf(KFE.s==0);
    stats.sav0 = sum(temp(:));

    lw_constrained_interp = griddedInterpolant(...
        grdKFE.b.vec, stats.cdf_b, 'pchip', 'nearest');

    w_constrained_interp = griddedInterpolant(...
        stats.values_cdf_wealth,...
        stats.cdf_wealth, 'pchip', 'nearest');

    stats.constrained_liq = lw_constrained_interp(p.epsilon_HtM);
    stats.constrained = w_constrained_interp(p.epsilon_HtM);

    %% --------------------------------------------------------------------
    % CONSTRAINED BY OWN QUARTERLY INCOME
    % ---------------------------------------------------------------------
    % ratio of total wealth to own quarterly income
    wi_ratio = (grdKFE.a.matrix + grdKFE.b.matrix) ./ income.y.matrixKFE;

    % ratio of liquid wealth to own quarterly income
    lwi_ratio = grdKFE.b.matrix ./ income.y.matrixKFE;

    % fraction with total wealth < own quarterly income / 6
    stats.HtM_one_sixth_Q_twealth = find_constrained(wi_ratio,stats.pmf, 1/6);
    % fraction with total wealth < own quarterly income / 12
    stats.HtM_one_twelfth_Q_twealth = find_constrained(wi_ratio,stats.pmf, 1/12);

    % fraction with liquid wealth < own quarterly income / 12
    stats.HtM_one_sixth_Q_lwealth = find_constrained(lwi_ratio,stats.pmf, 1/6);
    % fraction with liquid wealth < own quarterly income / 12
    stats.HtM_one_twelfth_Q_lwealth = find_constrained(lwi_ratio,stats.pmf, 1/12);

    stats.ratio_WHtM_HtM_sixth = 1 ...
    	- stats.HtM_one_sixth_Q_twealth / stats.HtM_one_sixth_Q_lwealth;
    stats.ratio_WHtM_HtM_twelfth = 1 ...
    	-  stats.HtM_one_twelfth_Q_twealth / stats.HtM_one_twelfth_Q_lwealth;
    
    %% --------------------------------------------------------------------
    % GINI COEFFICIENTS
    % ---------------------------------------------------------------------
    stats.wgini = direct_gini(stats.wealth_support, stats.pmf_wealth_support);

    %% --------------------------------------------------------------------
    % OUTPUT FOR HISTOGRAMS
    % ---------------------------------------------------------------------
    dstr = stats.pmf(:);
    [stats.b_hist.bins,stats.b_hist.values] = create_bins(1,grdKFE.b.matrix,dstr);
    [stats.a_hist.bins,stats.a_hist.values] = create_bins(1,grdKFE.a.matrix,dstr);
    [stats.w_hist.bins,stats.w_hist.values] = create_bins(1,wealth_mat,dstr);

    %% --------------------------------------------------------------------
    % ADJUSTMENT COSTS
    % ---------------------------------------------------------------------
    stats.adjcosts = struct();

    if ~p.OneAsset
        % derivative of adjustment costs?
        stats.adjcosts.chivar = p.chi1^(-p.chi2) / (1+p.chi2);

        % adjustment cost paid
        chii = HACTLib.aux.AdjustmentCost.cost(KFE.d(:),...
	        	grdKFE.a.matrix(:), p);
        
        % mean abs(d)/a
        d_div_a = abs(KFE.d(:) ./ max(grdKFE.a.matrix(:),p.a_lb));
        stats.adjcosts.mean_d_div_a = d_div_a' * stats.pmf(:);

        % median abs(d)/a
        stats.adjcosts.median_d_div_a = compute_pct(d_div_a, stats.pmf, 0.5);

        % condition on |d| > 0
        dvec = KFE.d(:);
        dnon0 = abs(dvec)>0;

        if sum(dnon0) > 0
        	% mean chi/abs(d) for |d| > 0
	        

	        chii_div_d = chii(dnon0) ./ abs(dvec(dnon0));
	        pmf_valid = stats.pmf(:);
	        pmf_valid = pmf_valid(dnon0) ./ sum(pmf_valid(dnon0));
	        stats.adjcosts.mean_chi_div_d = dot(chii_div_d, pmf_valid);

	        % median chi/abs(d) for |d| > 0
	        stats.adjcosts.median_chi_div_d = compute_pct(...
	        	chii_div_d, pmf_valid, 0.5);
	    else
	    	stats.adjcosts.mean_chi_div_d = NaN;
	    	stats.adjcosts.median_chi_div_d = NaN;
    	end

        % mean chi
        stats.adjcosts.mean_chi = dot(chii, stats.pmf(:));
        
        % fraction with d = 0
        stats.adjcosts.d0 = dot(KFE.d(:)==0,  stats.pmf(:));
    else
        stats.adjcosts.chivar = NaN;
        stats.adjcosts.mean_d_div_a = NaN;
        stats.adjcosts.median_d_div_a = NaN;
        stats.adjcosts.mean_chi_div_d = NaN;
        stats.adjcosts.median_chi_div_d = NaN;
        stats.adjcosts.mean_chi = NaN;
        stats.adjcosts.d0 = NaN;
    end

    %% --------------------------------------------------------------------
    % OTHER
    % ---------------------------------------------------------------------
    stats.meanc = dot(KFE.c(:), stats.pmf(:));
    
    % households that choose special case: to consume everything 
    % that's withdrawn
    stats.bmin_apos_consume_withdrawals = dot(...
    	KFE.bmin_consume_withdrawals(:), stats.pmf(:));
    
    bmin = grdKFE.b.matrix(:) == p.bmin;
    apos = grdKFE.a.matrix(:) > 0;
    bdot0 = abs(KFE.bdot(:)) < 1e-6;
    bdotpos = KFE.bdot(:) > 1e-6;
    adotpos = KFE.adot(:) > 0;
    adotneg = KFE.adot(:) < 0;
    dpos = KFE.d(:) > 0;
    dneg = KFE.d(:) < 0;
    d0 = KFE.d(:) == 0;
    
    % different possibilities for b = bmin, a > 0
    stats.bmin_apos = dot((bmin & apos), stats.pmf(:));
    
    stats.bmin_apos_adotpos = (bmin & apos & adotpos)' * stats.pmf(:);
    stats.bmin_apos_adotneg = (bmin & apos & adotneg)' * stats.pmf(:);
    
    stats.bmin_apos_bdotpos = (bmin & apos & bdotpos)' * stats.pmf(:);
    stats.bmin_apos_bdot0 = (bmin & apos & bdot0)' * stats.pmf(:);
    stats.bmin_apos_bdot0_dpos = (bmin & apos & bdot0 & dpos)' * stats.pmf(:);
    stats.bmin_apos_bdot0_dneg = (bmin & apos & bdot0 & dneg)' * stats.pmf(:);
    stats.bmin_apos_bdot0_d0 = (bmin & apos & bdot0 & d0)' * stats.pmf(:);
    
    stats.bmin_apos_bdot0_adotpos = (bmin & apos & bdot0 & adotpos)' * stats.pmf(:);
    stats.bmin_apos_bdot0_adotneg = (bmin & apos & bdot0 & adotneg)' * stats.pmf(:);
    
end


%% ------------------------------------------------------------------------
% FUNCTIONS
% -------------------------------------------------------------------------

%% Histogram helper function
function [bins,values] = create_bins(binwidth,vals,pdf1)
    sorts = sortrows([vals(:) pdf1(:)]);
    vals_sort = sorts(:,1);
    pmf = sorts(:,2);

    bins = min(vals(:)):binwidth:max(vals(:));
    values = zeros(size(bins));
    ibin = 1;
    for bn = bins
        if ibin == 1
            idx = vals_sort < bn + binwidth;
        elseif ibin < numel(bins)
            idx = (vals_sort>=bn) & (vals_sort<bn+binwidth);
        else
            idx = (vals_sort>=bn);
        end
        values(ibin) = sum(pmf(idx));
        ibin = ibin + 1;
    end
end