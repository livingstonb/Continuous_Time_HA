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
    % pmf(b)
    tmp = reshape(stats.pmf, p.nb_KFE, []);
    stats.pmf_b = sum(tmp, 2);

    b_support = stats.pmf_b > 1e-9;
    tmp = stats.pmf_b(b_support);
    stats.pmf_b_support = tmp / sum(tmp);
    stats.b_support = grdKFE.b.vec(b_support);

    tmp = cumsum(stats.pmf_b);
    stats.cdf_b_support = tmp(b_support);
    
    % pmf(a)
    tmp = reshape(stats.pmf, p.nb_KFE, p.na_KFE, []);
    stats.pmf_a = sum(sum(tmp, 3), 1);

    a_support = stats.pmf_a > 1e-9;
    tmp = stats.pmf_a(a_support);
    stats.pmf_a_support = tmp / sum(tmp);
    stats.a_support = grdKFE.a.vec(a_support);

    tmp = cumsum(stats.pmf_a);
    stats.cdf_a_support = tmp(a_support);
    
    % pmf(wealth)
    tmp = reshape(stats.pmf, p.nb_KFE*p.na_KFE, []);
    stats.pmf_wealth = sum(tmp, 2);
    wealth_support = stats.pmf_wealth > 1e-9;
    wealth_values = reshape(wealth_mat(:,:,1,1), [], 1);
    wealth_values_support = wealth_values(wealth_support);

    tmp = sort(wealth_values_support);
    unique_wealth_values = unique(tmp);
    n_wealth_support = numel(unique_wealth_values);

    stats.pmf_wealth_support = zeros(n_wealth_support, 1);
    for ii = 1:n_wealth_support
        wmask = (wealth_values == unique_wealth_values(ii));
        pmf_wealth_filtered = stats.pmf_wealth(wmask);
        stats.pmf_wealth_support(ii) = sum(pmf_wealth_filtered);
    end
    stats.pmf_wealth_support = stats.pmf_wealth_support...
        / sum(stats.pmf_wealth_support);
    stats.wealth_support = unique_wealth_values;

    stats.cdf_wealth_support = cumsum(stats.pmf_wealth_support);
  
    % Interpolants for percentiles
    lwinterp = griddedInterpolant(...
        stats.cdf_b_support,...
        stats.b_support,...
        'linear');

    if p.OneAsset == 0
        iwinterp = griddedInterpolant(...
            stats.cdf_a_support,...
            stats.a_support,...
            'linear');
    else
        iwinterp = @(a) NaN;
    end

    winterp = griddedInterpolant(...
        stats.cdf_wealth_support,...
        stats.wealth_support,...
        'linear');

    n_pct = numel(p.wpercentiles);
    stats.lwpercentile = zeros(1, n_pct);
    stats.iwpercentile = zeros(1, n_pct);
    stats.wpercentile = zeros(1, n_pct);
    for ii = 1:numel(p.wpercentiles)
        pct_val = p.wpercentiles(ii) / 100;
        stats.lwpercentile(ii) = lwinterp(pct_val);
        stats.iwpercentile(ii) = iwinterp(pct_val);
        stats.wpercentile(ii) = winterp(pct_val);
    end

    stats.median_liqw = lwinterp(0.5);
    stats.median_illliqw = iwinterp(0.5);
    stats.median_totw = winterp(0.5);

    % Top shares
    % Amount of total assets that reside in each pt on sorted asset space
    totassets = stats.pmf_wealth_support .* stats.wealth_support;
    % Fraction of total assets in each pt on asset space
    cumassets = cumsum(totassets) / stats.totw;
    cumwealthshare_interp = griddedInterpolant(...
        stats.cdf_wealth_support,...
        cumassets,...
        'linear');
    stats.top10share = 1 - cumwealthshare_interp(0.9);
    stats.top1share = 1 - cumwealthshare_interp(0.99);

    %% --------------------------------------------------------------------
    % CONSTRAINED HOUSEHOLDS
    % ---------------------------------------------------------------------
    n_HtM = numel(p.epsilon_HtM);
    stats.constrained = zeros(1, n_HtM);
    stats.constrained_liq = zeros(1, n_HtM);

    % total wealth = 0
    temp = stats.pmf(wealth_mat==0);
    stats.constrained(1) = sum(temp(:));

    % liquid wealth = 0
    temp = stats.pmf(grdKFE.b.matrix==0);
    stats.constrained_liq(1) = sum(temp(:));

    % saving = 0
    temp = stats.pmf(KFE.s==0);
    stats.sav0 = sum(temp(:));

    lw_constrained_interp = griddedInterpolant(...
        stats.b_support,...
        stats.cdf_b_support,...
        'linear');

    w_constrained_interp = griddedInterpolant(...
        stats.wealth_support,...
        stats.cdf_wealth_support,...
        'linear');

    for ii = 1:n_HtM
        threshold = p.epsilon_HtM(ii);

        if threshold == 0
            temp = stats.pmf(wealth_mat==0);
            stats.constrained(ii) = sum(temp(:));

            temp = stats.pmf(grdKFE.b.matrix==0);
            stats.constrained_liq(ii) = sum(temp(:));
        else
            stats.constrained(ii) = w_constrained_interp(threshold);
            stats.constrained_liq(ii) = lw_constrained_interp(threshold);
        end
    end

    %% --------------------------------------------------------------------
    % CONSTRAINED BY OWN QUARTERLY INCOME
    % ---------------------------------------------------------------------
    % ratio of total wealth to own quarterly income
    wi_ratio = (grdKFE.a.matrix + grdKFE.b.matrix) ./ income.y.matrixKFE;

    % ratio of liquid wealth to own quarterly income
    lwi_ratio = grdKFE.b.matrix ./ income.y.matrixKFE;

    % fraction with total wealth < own quarterly income / 6
    stats.HtM_one_sixth_Q_twealth = find_constrained(wi_ratio,stats.pmf,1/6);
    % fraction with total wealth < own quarterly income / 12
    stats.HtM_one_twelfth_Q_twealth = find_constrained(wi_ratio,stats.pmf,1/12);

    % fraction with liquid wealth < own quarterly income / 12
    stats.HtM_one_sixth_Q_lwealth = find_constrained(lwi_ratio,stats.pmf,1/6);
    % fraction with liquid wealth < own quarterly income / 12
    stats.HtM_one_twelfth_Q_lwealth = find_constrained(lwi_ratio,stats.pmf,1/12);

    stats.ratio_WHtM_HtM_sixth = 1 -  stats.HtM_one_sixth_Q_twealth / stats.HtM_one_sixth_Q_lwealth;
    stats.ratio_WHtM_HtM_twelfth = 1 -  stats.HtM_one_twelfth_Q_twealth / stats.HtM_one_twelfth_Q_lwealth;
    
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

    if p.OneAsset == 0
        % derivative of adjustment costs?
        stats.adjcosts.chivar = p.chi1^(-p.chi2) / (1+p.chi2);
        
        % mean abs(d)/a
        d_div_a = abs(KFE.d(:) ./ max(grdKFE.a.matrix(:),p.a_lb));
        stats.adjcosts.mean_d_div_a = d_div_a' * stats.pmf(:);

        % median abs(d)/a
        stats.adjcosts.median_d_div_a = compute_pct(d_div_a,stats.pmf,0.5);

        % mean chi/abs(d) for |d| > 0
        chii = HACTLib.aux.AdjustmentCost.cost(KFE.d(:),grdKFE.a.matrix(:),p);
        chii_div_d = chii ./ abs(KFE.d(:));
        pmf_valid = stats.pmf(:);
        pmf_valid = pmf_valid(abs(KFE.d(:))>0) ./ sum(pmf_valid(abs(KFE.d(:))>0));
        stats.adjcosts.mean_chi_div_d = chii_div_d(abs(KFE.d(:))>0)' * pmf_valid;

        % median chi/abs(d) for d > 0 
        chid_sort = sortrows([chii_div_d(abs(KFE.d(:))>0) pmf_valid]);
        chid_values_sort = chid_sort(:,1);
        chid_cdf = cumsum(chid_sort(:,2));
        [chid_cdf_u,uind] = unique(chid_cdf,'last');
        chid_values_u = chid_values_sort(uind);
        try
            chid_interp = griddedInterpolant(chid_cdf_u,chid_values_u,'linear');
            stats.adjcosts.median_chi_div_d = chid_interp(0.5);
        catch
            stats.adjcosts.median_chi_div_d = NaN;
        end

        % mean chi
        stats.adjcosts.mean_chi = chii' * stats.pmf(:);
        
        % fraction with d = 0
        stats.adjcosts.d0 = (KFE.d(:)==0)' * stats.pmf(:);
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
    stats.meanc = KFE.c(:)' * stats.pmf(:);
    
    % households that choose special case: to consume everything 
    % that's withdrawn
    stats.bmin_apos_consume_withdrawals = KFE.bmin_consume_withdrawals(:)' * stats.pmf(:);
    
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
    stats.bmin_apos = (bmin & apos)' * stats.pmf(:);
    
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