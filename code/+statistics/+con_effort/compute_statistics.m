function stats = compute_statistics(p,income,grdKFE,KFE)
    %% --------------------------------------------------------------------
    % Consumption by hand-to-mouth status
    %
    % KFE.c is just the consumption grid in smooth adjustment case.
    % With instantaneous switching, KFE.c is adjusted based on
    % which state is switched to.
    % ---------------------------------------------------------------------

	stats.pmf = KFE.g .* grdKFE.trapezoidal.matrix;
    stats.rho = p.rho;

    %% --------------------------------------------------------------------
    % WEALTH LEVELS
    % ---------------------------------------------------------------------
    stats.wealth = stats.pmf(:)' * grdKFE.b.matrix(:);
    
    %% --------------------------------------------------------------------
    % WEALTH PERCENTILES
    % ---------------------------------------------------------------------

    % total wealth
    stats.wpercentile = aux.compute_pct(...
    	grdKFE.b.matrix,stats.pmf,p.wpercentiles/100);
    
    % top wealth shares
    wsort = sortrows([grdKFE.b.matrix(:) stats.pmf(:)]);
    cdf_wealth = cumsum(wsort(:,2));
    [cdf_wealth_u,uind] = unique(cdf_wealth,'last');
    % amount of total assets that reside in each point on asset space
    totassets = wsort(:,2) .* wsort(:,1);
    % cumulative fraction of total assets at each point in asset space
    cumassets = cumsum(totassets) / stats.wealth;
    % create interpolant from wealth percentile to cumulative wealth share
    cumwealthshare = griddedInterpolant(cdf_wealth_u,cumassets(uind),'linear');
    stats.top10share = 1 - cumwealthshare(0.9);
    stats.top1share = 1 - cumwealthshare(0.99);

    %% --------------------------------------------------------------------
    % CONSTRAINED HOUSEHOLDS
    % ---------------------------------------------------------------------
    % wealth == 0
    temp = stats.pmf(grdKFE.b.matrix==0);
    stats.constrained(1) = sum(temp(:));

    % wealth < epsilon
    stats.constrained(2:numel(p.epsilon_HtM)) = ...
    	aux.find_constrained(grdKFE.b.matrix,stats.pmf,p.epsilon_HtM(2:end));

    %% --------------------------------------------------------------------
    % CONSTRAINED BY OWN QUARTERLY INCOME
    % ---------------------------------------------------------------------
    % construct cdf of total wealth / own quarterly income
    wi_ratio = grdKFE.b.matrix ./ income.y.matrixKFE;

    % fraction with total wealth < own quarterly income / 6
    stats.HtM_one_sixth_Q_wealth = find_constrained(wi_ratio,stats.pmf,1/6);
    % fraction with total wealth < own quarterly income / 12
    stats.HtM_one_twelfth_Q_wealth = find_constrained(wi_ratio,stats.pmf,1/12);

    %% --------------------------------------------------------------------
    % CONSUMPTION
    % ---------------------------------------------------------------------
    % |h|
    stats.mean_absh_total = abs(KFE.h(:))' * stats.pmf(:);
    hnonzero = abs(KFE.h)>1e-8;
    hnonzeromass = sum(reshape(stats.pmf(hnonzero),[],1));
    stats.mean_absh_ofchangers = abs(KFE.h(hnonzero))'...
                                    * stats.pmf(hnonzero) / hnonzeromass;
                                
    % sign of h
    stats.h0 = (abs(KFE.h(:))<1e-8)' * stats.pmf(:);
    stats.hpos = (KFE.h(:)>0)' * stats.pmf(:);
    stats.hneg = (KFE.h(:)<0)' * stats.pmf(:);
    
    % consumption distribution
    stats.c10 = aux.compute_pct(KFE.c,stats.pmf,0.1);
    stats.c95 = aux.compute_pct(KFE.c,stats.pmf,0.95);
    stats.c99 = aux.compute_pct(KFE.c,stats.pmf,0.99);
    stats.c999 = aux.compute_pct(KFE.c,stats.pmf,0.999);

    %% --------------------------------------------------------------------
    % FRACTION OF HHs WITH a < 0
    % ---------------------------------------------------------------------
    stats.anegative = (grdKFE.b.matrix(:)<0)' * stats.pmf(:);

    %% --------------------------------------------------------------------
    % GINI COEFFICIENTS
    % ---------------------------------------------------------------------
    stats.wgini = direct_gini(grdKFE.b.matrix,stats.pmf);

    %% --------------------------------------------------------------------
    % OUTPUT FOR HISTOGRAMS (NOT USED)
    % ---------------------------------------------------------------------
    [stats.a_hist.bins,stats.a_hist.values] = create_bins(1,grdKFE.b.matrix,stats.pmf(:));
    
    %% --------------------------------------------------------------------
    % OTHER
    % ---------------------------------------------------------------------
    stats.beta_annualized = exp(-4 * p.rho);
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