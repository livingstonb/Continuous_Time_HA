clearvars -except stats p

%% This script is used to combine one or more variablesX.mat files. Produces a table.

%% Set FROM_MATFILE = false if running right after model, true if running from .mat file
FROM_MATFILE = true;

%% Select directories
matdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output/';
% matdir = '/home/brian/Documents/GitHub/Continuous_Time_HA/output/';
% matdir = '/media/hdd/Other/midway2_output/continuous_time';
% matdir = '/Users/brianlivingston/Documents/midway2_output/';
% matdir = '/Users/Brian-laptop/Documents/GitHub/Continuous_Time_HA/output/';

codedir = '/home/livingstonb/GitHub/Continuous_Time_HA/code/';
% codedir = '/home/brian/Documents/GitHub/Continuous_Time_HA/code/';

% matdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output/';
% codedir = '/home/livingstonb/GitHub/Continuous_Time_HA/';
xlxpath = '/home/livingstonb/GitHub/Continuous_Time_HA/output/';
xlxpath1 = [xlxpath 'detailedResults.xlsx'];
xlxpath2 = [xlxpath 'decomposition.xlsx'];

addpath(codedir);

%% read baseline one-asset for decomposition
fpath = [matdir,'output_oneasset.mat'];
if exist(fpath,'file')
    oneasset = load(fpath);
else
    disp('no one asset baseline found')
    oneasset = [];
end

%% Read .mat files into a cell array
if FROM_MATFILE
    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [matdir,'output_',runstr,'.mat'];
        if exist(fpath,'file')
            ind = ind+1;

            % to save memory
            s(ind) = load(fpath);

            if ind > 1
                s(ind).grd = [];
                s(ind).KFE = [];
            end

            % perform Empc1 - Empc0 decomposition
            decomp_base(ind) = statistics.decomp_baseline(s(1),s(ind));    

            % perform decomp wrt one-asset model
            decomp_oneasset(ind) = statistics.decomp_twoasset_oneasset(oneasset,s(ind));
        else
            continue
        end
    end
else
    s.p = p;
    s.stats = stats;
    
    skip = true;
    decomp_base = statistics.decomp_baseline(s(1),s(1));
    decomp_oneasset = statistics.decomp_twoasset_oneasset(oneasset,s(1));
end

n = numel(s);
nans = num2cell(NaN(1,n));

import HACTLib.aux.get_all_values

detailedResults = [{'Name'}, get_all_values(s,'p',1,'name')
            {'chi0'}, get_all_values(s,'p',1,'chi0')
            {'chi1'}, get_all_values(s,'p',1,'chi1')
            {'chi2'}, get_all_values(s,'p',1,'chi2')
            {'chi1^(-chi2)/(1+chi2)'}, get_all_values(s,'stats',1,'adjcosts',1,'chivar')
            {'Rho'}, get_all_values(s,'p',1,'rho')
            {'Beta (Annualized)'}, get_all_values(s,'stats',1,'beta_annualized')
            {'r_a'}, get_all_values(s,'p',1,'r_a')
            {'Illiq Assets'}, get_all_values(s,'stats',1,'illiqw')
            {'Liq Assets'}, get_all_values(s,'stats',1,'liqw')
            {'Total Assets'}, get_all_values(s,'stats',1,'totw')
            {'Wealth Top 10% Share'}, get_all_values(s,'stats',1,'top10share')
            {'Wealth Top 1% Share'}, get_all_values(s,'stats',1,'top1share')
            {'Gini (Total Assets)'}, get_all_values(s,'stats',1,'wgini')
            {'s = 0'}, get_all_values(s,'stats',1,'sav0')
            {'____HtM STATISTICS'}, nans
            {'Wealth HtM / Total HtM (at 1/6 qincome)'}, get_all_values(s,'stats',1,'ratio_WHtM_HtM_sixth');
            {'Wealth HtM / Total HtM (at 1/12 qincome)'}, get_all_values(s,'stats',1,'ratio_WHtM_HtM_twelfth');
            {'____CONSTRAINED IN LIQUID WEALTH'}, nans
            {'L Wealth == 0'}, get_all_values(s,'stats',1,'constrained_liq',1)
            {'L Wealth <= 0.5% Mean Annual Income'}, get_all_values(s,'stats',1,'constrained_liq',2)
            {'L Wealth <= 1% Mean Annual Income'}, get_all_values(s,'stats',1,'constrained_liq',3)
            {'L Wealth <= 2% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained_liq',4)
            {'L Wealth <= 5% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained_liq',5)
            {'L Wealth <= 10% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained_liq',6)
            {'L Wealth <= 15% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained_liq',7)
            {'L Wealth <= (1/6) * Own Quarterly Income'}, get_all_values(s,'stats',1,'HtM_one_sixth_Q_lwealth')
            {'L Wealth <= (1/12) * Own Quarterly Income'}, get_all_values(s,'stats',1,'HtM_one_twelfth_Q_lwealth')
            {'____CONSTRAINED IN TOTAL WEALTH'}, nans
            {'Wealth == 0'}, get_all_values(s,'stats',1,'constrained',1)
            {'Wealth <= 0.5% Mean Annual Income'}, get_all_values(s,'stats',1,'constrained',2)
            {'Wealth <= 1% Mean Annual Income'}, get_all_values(s,'stats',1,'constrained',3)
            {'Wealth <= 2% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained',4)
            {'Wealth <= 5% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained',5)
            {'Wealth <= 10% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained',6)
            {'Wealth <= 15% * Mean Annual Income'}, get_all_values(s,'stats',1,'constrained',7)
            {'Wealth <= (1/6) * Own Quarterly Income'}, get_all_values(s,'stats',1,'HtM_one_sixth_Q_twealth')
            {'Wealth <= (1/12) * Own Quarterly Income'}, get_all_values(s,'stats',1,'HtM_one_twelfth_Q_twealth')
            {'____ILLIQUID ASSETS ADJUSTMENT'}, nans
			{'Mean |d|/max(a,a_lb)'}, get_all_values(s,'stats',1,'adjcosts',1,'mean_d_div_a')
            {'Median |d|/max(a,a_lb)'}, get_all_values(s,'stats',1,'adjcosts',1,'median_d_div_a')
            {'Mean chi/|d|, given |d|>0'}, get_all_values(s,'stats',1,'adjcosts',1,'mean_chi_div_d')
            {'Median chi/|d|, given |d|>0'}, get_all_values(s,'stats',1,'adjcosts',1,'median_chi_div_d')
            {'Mean chi'}, get_all_values(s,'stats',1,'adjcosts',1,'mean_chi')
            {'Fraction with d == 0'}, get_all_values(s,'stats',1,'adjcosts',1,'d0')
            {'____WEALTH PERCENTILES'}, nans
            {'Liquid wealth, 10th'}, get_all_values(s,'stats',1,'lwpercentile',1)
            {'Liquid wealth, 25th'}, get_all_values(s,'stats',1,'lwpercentile',2)
            {'Liquid wealth, 50th'}, get_all_values(s,'stats',1,'lwpercentile',3)
            {'Liquid wealth, 90th'}, get_all_values(s,'stats',1,'lwpercentile',4)
            {'Liquid wealth, 99th'}, get_all_values(s,'stats',1,'lwpercentile',5)
            {'Liquid wealth, 99.9th'}, get_all_values(s,'stats',1,'lwpercentile',6)
            {'Total wealth, 10th'}, get_all_values(s,'stats',1,'wpercentile',1)
            {'Total wealth, 25th'}, get_all_values(s,'stats',1,'wpercentile',2)
            {'Total wealth, 50th'}, get_all_values(s,'stats',1,'wpercentile',3)
            {'Total wealth, 90th'}, get_all_values(s,'stats',1,'wpercentile',4)
            {'Total wealth, 99th'}, get_all_values(s,'stats',1,'wpercentile',5)
            {'Total wealth, 99.9th'}, get_all_values(s,'stats',1,'wpercentile',6)
            {'____AVG MPC OUT OF -1e-5 * MEAN ANNUAL INCOME'}, nans
            {'QUARTER 1 MPC, shock = -1e-5'}, get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',1)
            {'QUARTER 2 MPC, shock = -1e-5'}, get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',2)
            {'QUARTER 3 MPC shock = -1e-5'}, get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',3)
            {'QUARTER 4 MPC shock = -1e-5'}, get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',4)
            {'ANNUAL MPC, shock = -1e-5'}, get_all_values(s,'stats',1,'mpcs',1,'avg_0_annual')
            {'____AVG MPC OUT OF -0.01 * MEAN ANNUAL INCOME'}, nans
            {'QUARTER 1 MPC, shock = -0.01'}, get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',1)
            {'QUARTER 2 MPC, shock = -0.01'}, get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',2)
            {'QUARTER 3 MPC shock = -0.01'}, get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',3)
            {'QUARTER 4 MPC shock = -0.01'}, get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',4)
            {'ANNUAL MPC, shock = -0.01'}, get_all_values(s,'stats',1,'mpcs',2,'avg_0_annual')
            {'____AVG MPC OUT OF -0.1 * MEAN ANNUAL INCOME'}, nans
            {'QUARTER 1 MPC, shock = -0.1'}, get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',1)
            {'QUARTER 2 MPC, shock = -0.1'}, get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',2)
            {'QUARTER 3 MPC shock = -0.1'}, get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',3)
            {'QUARTER 4 MPC shock = -0.1'}, get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',4)
            {'ANNUAL MPC, shock = -0.1'}, get_all_values(s,'stats',1,'mpcs',3,'avg_0_annual')
            {'____AVG MPC OUT OF 1e-5 * MEAN ANNUAL INCOME'}, nans
            {'QUARTER 1 MPC, shock = 1e-5'}, get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',1)
            {'QUARTER 2 MPC, shock = 1e-5'}, get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',2)
            {'QUARTER 3 MPC shock = 1-5'}, get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',3)
            {'QUARTER 4 MPC shock = 1e-5'}, get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',4)
            {'ANNUAL MPC, shock = 1e-5'}, get_all_values(s,'stats',1,'mpcs',4,'avg_0_annual')
            {'____AVG MPC OUT OF 0.01 * MEAN ANNUAL INCOME'}, nans
            {'QUARTER 1 MPC, shock = 0.01'}, get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',1)
            {'QUARTER 2 MPC, shock = 0.01'}, get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',2)
            {'QUARTER 3 MPC shock = 0.01'}, get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',3)
            {'QUARTER 4 MPC shock = 0.01'}, get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',4)
            {'ANNUAL MPC, shock = 0.01'}, get_all_values(s,'stats',1,'mpcs',5,'avg_0_annual')
            {'____AVG MPC OUT OF 0.1 * MEAN ANNUAL INCOME'}, nans
            {'QUARTER 1 MPC, shock = 0.1'}, get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',1)
            {'QUARTER 2 MPC, shock = 0.1'}, get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',2)
            {'QUARTER 3 MPC shock = 0.1'}, get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',3)
            {'QUARTER 4 MPC shock = 0.1'}, get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',4)
            {'ANNUAL MPC, shock = 0.1'}, get_all_values(s,'stats',1,'mpcs',6,'avg_0_annual')
            {'____AVG QUARTERLY MPC OUT OF NEWS'}, nans
            {'QUARTER 1 MPC, shock = -0.01 next quarter'}, get_all_values(s,'stats',1,'mpcs',2,'avg_1_quarterly')
            {'ANNUAL MPC, shock = -0.01 next year'}, get_all_values(s,'stats',1,'mpcs',2,'avg_4_annual',1)
            {'QUARTER 1 MPC, shock = 0.01 next quarter'}, get_all_values(s,'stats',1,'mpcs',5,'avg_1_quarterly')
            {'QUARTER 1 MPC, shock = 0.1 next quarter'}, get_all_values(s,'stats',1,'mpcs',6,'avg_1_quarterly')
            {'QUARTER 1 MPC, shock = 0.01 next year'}, get_all_values(s,'stats',1,'mpcs',5,'avg_4_quarterly',1)
            {'QUARTER 1 MPC, shock = 0.1 next year'}, get_all_values(s,'stats',1,'mpcs',6,'avg_4_quarterly',1)
            {'ANNUAL MPC, shock = 0.01 next year'}, get_all_values(s,'stats',1,'mpcs',5,'avg_4_annual')
            {'ANNUAL MPC, shock = 0.1 next year'}, get_all_values(s,'stats',1,'mpcs',6,'avg_4_annual')
            {'____SOME SIMULATED MEAN MPCs'}, nans
            {'QUARTER 1 SIM MPC, shock = 0.01'}, get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly',1)
            {'QUARTER 1 SIM MPC, shock = 0.1'}, get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',1)
            {'___DECOMP OF EM1 AROUND 0'}, nans
            {'Decomp around 0, RA MPC'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term1',1)
            {'Decomp around 0, HtM Effect'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term2',1)
            {'Decomp around 0, Non-HtM, constraint'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term3',1)
            {'Decomp around 0, Non-HtM, inc risk'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term4',1)
            {'___DECOMP OF EM1 AROUND 0.01'}, nans
            {'Decomp around 0.01, RA MPC'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term1',2)
            {'Decomp around 0.01, HtM Effect'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term2',2)
            {'Decomp around 0.01, Non-HtM, constraint'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term3',2)
            {'decomp Around 0.01, Non-HtM, inc risk'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term4',2)
            {'___DECOMP OF EM1 AROUND 0.05'}, nans
            {'Decomp around 0.05, RA MPC'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term1',3)
            {'Decomp around 0.05, HtM Effect'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term2',3)
            {'Decomp around 0.05, Non-HtM, constraint'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term3',3)
            {'Decomp around 0.05, Non-HtM, inc risk'}, get_all_values(s,'stats',1,'decomp_norisk',1,'term4',3)
            {'___DECOMP OF EM1-EM0 (MPC OUT OF 0.01)'}, nans
            {'Em1 - Em0'}, get_all_values(decomp_base,'Em1_less_Em0')
            {'Em1 - Em0, Effect of MPC fcn'}, get_all_values(decomp_base,'term1')
            {'Em1 - Em0, Effect of distr'}, get_all_values(decomp_base,'term2')
            {'Em1 - Em0, Interaction'}, get_all_values(decomp_base,'term3')
            {'Distr effect around 0, PHtM only'}, get_all_values(decomp_base,'term2a',1)
            {'Distr effect around 0, WHtM only'}, get_all_values(decomp_base,'term2b',1)
            {'Distr effect around 0, non-HtM only'}, get_all_values(decomp_base,'term2c',1)
            {'Distr effect around 0.01, PHtM only'}, get_all_values(decomp_base,'term2a',2)
            {'Distr effect around 0.01, WHtM only'}, get_all_values(decomp_base,'term2b',2)
            {'Distr effect around 0.01, non-HtM only'}, get_all_values(decomp_base,'term2c',2)
            {'Distr effect around 0.05, PHtM only'}, get_all_values(decomp_base,'term2a',3)
            {'Distr effect around 0.05, WHtM only'}, get_all_values(decomp_base,'term2b',3)
            {'Distr effect around 0.05, non-HtM only'}, get_all_values(decomp_base,'term2c',3)
            {'___DECOMP WRT ONE-ASSET MODEL'}, nans
            {'E[mpc_twoasset] - E[mpc_oneasset]'}, get_all_values(decomp_oneasset,'Em1_minus_Em0')
            {'Effect of MPC fn'}, get_all_values(decomp_oneasset,'term1_mpcfn')
            {'Effect of MPC fn, PHtM term around 0'}, get_all_values(decomp_oneasset,'term1a',1)
            {'Effect of MPC fn, WHtM term around 0'}, get_all_values(decomp_oneasset,'term1b',1)
            {'Effect of MPC fn, NHtM term around 0'}, get_all_values(decomp_oneasset,'term1c',1)
            {'Effect of MPC fn, PHtM term around 0.01'}, get_all_values(decomp_oneasset,'term1a',2)
            {'Effect of MPC fn, WHtM term around 0.01'}, get_all_values(decomp_oneasset,'term1b',2)
            {'Effect of MPC fn, NHtM term around 0.01'}, get_all_values(decomp_oneasset,'term1c',2)
            {'Effect of distr of net worth'}, get_all_values(decomp_oneasset,'term2_networth')
            {'Effect of distr of net worth, PHtM term around 0'}, get_all_values(decomp_oneasset,'term2a',1)
            {'Effect of distr of net worth, n>0 term around 0'}, get_all_values(decomp_oneasset,'term2b',1)
            {'Effect of distr of net worth, PHtM term around 0.01'}, get_all_values(decomp_oneasset,'term2a',2)
            {'Effect of distr of net worth, n>0 term around 0.01'}, get_all_values(decomp_oneasset,'term2b',2)
            {'Interaction'}, get_all_values(decomp_oneasset,'term3_interaction')
            {'___DECOMP OF EM1-RA MPC'}, nans
            {'Em1 - EmRA'}, get_all_values(s,'stats',1,'decompRA',1,'Em1_less_mRA')
            {'Em1 - EmRA, Effect of MPC fcn'}, get_all_values(s,'stats',1,'decompRA',1,'term1')
            {'Em1 - EmRA, Effect of distr'}, get_all_values(s,'stats',1,'decompRA',1,'term2')
            {'Em1 - EmRA. Interaction'}, get_all_values(s,'stats',1,'decompRA',1,'term3')];

Em1_less_Em0 = cell2mat(get_all_values(decomp_base,'Em1_less_Em0'));
decompTable = [	{'Name'}, get_all_values(s,'p',1,'name')
				{'Em1 - Em0'}, get_all_values(decomp_base,'Em1_less_Em0')
				{'Effect of MPC fcn'}, get_all_values(decomp_base,'term1')
				{'Effect of distr'}, get_all_values(decomp_base,'term2')
				{'Distr effect around 0, PHtM only'}, get_all_values(decomp_base,'term2a',1)
	            {'Distr effect around 0, WHtM only'}, get_all_values(decomp_base,'term2b',1)
	            {'Distr effect around 0, non-HtM only'}, get_all_values(decomp_base,'term2c',1)
	            {'Distr effect around 0.01, PHtM only'}, get_all_values(decomp_base,'term2a',2)
	            {'Distr effect around 0.01, WHtM only'}, get_all_values(decomp_base,'term2b',2)
	            {'Distr effect around 0.01, non-HtM only'}, get_all_values(decomp_base,'term2c',2)
            	{'Em1 - Em0, Interaction'}, get_all_values(decomp_base,'term3')
            	{'% Effect of MPC fcn'}, num2cell(cell2mat(get_all_values(decomp_base,'term1')) ./ Em1_less_Em0 * 100)
				{'% Effect of distr'}, num2cell(cell2mat(get_all_values(decomp_base,'term2')) ./ Em1_less_Em0 * 100)
				{'% Distr effect around 0, PHtM only'}, num2cell(cell2mat(get_all_values(decomp_base,'term2a',1)) ./ Em1_less_Em0 * 100)
	            {'% Distr effect around 0, WHtM only'}, num2cell(cell2mat(get_all_values(decomp_base,'term2b',1)) ./ Em1_less_Em0 * 100)
	            {'% Distr effect around 0, non-HtM only'}, num2cell(cell2mat(get_all_values(decomp_base,'term2c',1)) ./ Em1_less_Em0 * 100)
	            {'% Distr effect around 0.01, PHtM only'}, num2cell(cell2mat(get_all_values(decomp_base,'term2a',2)) ./ Em1_less_Em0 * 100)
	            {'% Distr effect around 0.01, WHtM only'}, num2cell(cell2mat(get_all_values(decomp_base,'term2b',2)) ./ Em1_less_Em0 * 100)
	            {'% Distr effect around 0.01, non-HtM only'}, num2cell(cell2mat(get_all_values(decomp_base,'term2c',2)) ./ Em1_less_Em0 * 100)
            	{'% Em1 - Em0, Interaction'}, num2cell(cell2mat(get_all_values(decomp_base,'term3')) ./ Em1_less_Em0 * 100)];

TdetailedResults = cell2table(detailedResults(:,2:n+1),'RowNames',detailedResults(:,1));
TdecompTable = cell2table(decompTable(:,2:n+1),'RowNames',decompTable(:,1));

if ~isempty(xlxpath)
	writetable(TdetailedResults,xlxpath1,'WriteRowNames',true)
	writetable(TdecompTable,xlxpath2,'WriteRowNames',true)
end
