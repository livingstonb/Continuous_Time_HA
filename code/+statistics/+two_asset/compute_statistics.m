function stats = statistics(p,income,grd,grdKFE,KFE)
	% computes various statistics after solving the model
	
    %% --------------------------------------------------------------------
    % Consumption by hand-to-mouth status
    % ---------------------------------------------------------------------
	na_KFE = p.na_KFE;
	nb_KFE = p.nb_KFE;
	ny = p.ny;
    nz = p.nz;

	stats.pmf = KFE.g .* grdKFE.trapezoidal.matrix;
    stats.pmf_norisk = sum(stats.pmf,3);
    stats.rho = p.rho;
    stats.beta = exp(-4*p.rho);
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

    stats.wpercentile = compute_pct(wealth_mat,stats.pmf,p.wpercentiles/100);
    stats.lwpercentile = compute_pct(grdKFE.b.matrix,stats.pmf,p.wpercentiles/100);
    stats.iwpercentile = compute_pct(grdKFE.a.matrix,stats.pmf,p.wpercentiles/100);

    % top wealth shares
    wsort = sortrows([wealth_mat(:) stats.pmf(:)]);
    cdf_wealth = cumsum(wsort(:,2));
    [cdf_wealth_u,uind] = unique(cdf_wealth,'last');
    % amount of total assets that reside in each point on asset space
    totassets = wsort(:,2) .* wsort(:,1);
    % cumulative fraction of total assets at each point in asset space
    cumassets = cumsum(totassets) / stats.totw;
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(cdf_wealth_u,cumassets(uind),'linear');
    stats.top10share = 1 - cumwealthshare(0.9);
    stats.top1share = 1 - cumwealthshare(0.99);

    %% --------------------------------------------------------------------
    % CONSTRAINED HOUSEHOLDS
    % ---------------------------------------------------------------------

    % total wealth = 0
    temp = stats.pmf(wealth_mat==0);
    stats.constrained(1) = sum(temp(:));

    % liquid wealth = 0
    temp = stats.pmf(grdKFE.b.matrix==0);
    stats.constrained_liq(1) = sum(temp(:));

    % saving = 0
    temp = stats.pmf(KFE.s==0);
    stats.sav0 = sum(temp(:));

    % wealth < epsilon
    stats.constrained(2:numel(p.epsilon_HtM)) = ...
    	find_constrained(wealth_mat,stats.pmf,p.epsilon_HtM(2:end));

    % liquid wealth < epsilon
    stats.constrained_liq(2:numel(p.epsilon_HtM)) = ...
    	find_constrained(grdKFE.b.matrix,stats.pmf,p.epsilon_HtM(2:end));

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
    % fraction with total wealth < own annual income * 0.15
    stats.HtM_015_ann_inc_twealth = find_constrained(wi_ratio,stats.pmf,0.15*4);

    % fraction with liquid wealth < own quarterly income / 12
    stats.HtM_one_sixth_Q_lwealth = find_constrained(lwi_ratio,stats.pmf,1/6);
    % fraction with liquid wealth < own quarterly income / 12
    stats.HtM_one_twelfth_Q_lwealth = find_constrained(lwi_ratio,stats.pmf,1/12);
    % fraction with liquid wealth < own annual income * 0.15
    stats.HtM_015_ann_inc_lwealth = find_constrained(lwi_ratio,stats.pmf,0.15*4);
    
    %% --------------------------------------------------------------------
    % GINI COEFFICIENTS
    % ---------------------------------------------------------------------
    stats.wgini = direct_gini(wealth_mat,stats.pmf);

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
        chii = aux.two_asset.adj_cost_fn(KFE.d(:),grdKFE.a.matrix(:),p);
        chii_div_d = chii ./ abs(KFE.d(:));
        ptmass_valid = stats.pmf(:);
        ptmass_valid = ptmass_valid(abs(KFE.d(:))>0) ./ sum(ptmass_valid(abs(KFE.d(:))>0));
        stats.adjcosts.mean_chi_div_d = chii_div_d(abs(KFE.d(:))>0)' * ptmass_valid;

        % median chi/abs(d) for d > 0 
        chid_sort = sortrows([chii_div_d(abs(KFE.d(:))>0) ptmass_valid]);
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
    stats.beta = exp(-4 * p.rho);
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

function gini = direct_gini(level,distr)
    % Sort distribution and levels by levels
    sort1 = sortrows([level(:),distr(:)]);
    level_sort = sort1(:,1);
    dist_sort  = sort1(:,2);
    S = [0;cumsum(dist_sort .* level_sort)];
    gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
end

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

function percentiles_values = compute_pct(values,pmf,percentiles)
	% finds percentiles, given a pmf over values

	% 'percentiles' input vector should be fractions, i.e.
	% to compute 99th percentile, use 0.99
	
	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[cdf_unique,unique_indices] = unique(cdf_sorted,'last');
	values_sorted_unique = values_sorted(unique_indices);

	cdf_interp = griddedInterpolant(cdf_unique,values_sorted_unique,'linear');
	percentiles_values = cdf_interp(percentiles);
end

function fraction_below = find_constrained(values,pmf,thresholds);
	% given values and a pmf, this functions finds the fraction
	% of households below given thresholds
	sorted_by_values = sortrows([values(:) pmf(:)]);
	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[values_unique,unique_indices] = unique(values_sorted,'last');
	cdf_unique = cdf_sorted(unique_indices);

	values_interp = griddedInterpolant(values_unique,cdf_unique,'linear');
	fraction_below = values_interp(thresholds);
end