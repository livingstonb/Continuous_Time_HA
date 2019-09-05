clear all

% This script is used to combine the variablesX.mat files produced when run
% in batch for the continuous time model. Produces a table with the variable
% name 'T'. For the decompositions, variables1.mat is used as the baseline.

% User must set matdir as location of the variablesX.mat files, where X is
% a number between 1 and 999.

%% Select directories
% matdir = '/home/brian/Documents/GitHub/Continuous_Two_Asset/Output';
% matdir = '/media/hdd/Other/midway2_output/continuous_time';
% matdir = '/Users/Brian-laptop/Documents/midway2_output/8_2_19/';
% matdir = '/Users/Brian-laptop/Documents/GitHub/Continuous_Two_Asset/Output/';

matdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output/two_asset/';
codedir = '/home/livingstonb/GitHub/Continuous_Time_HA/';
xlxpath = '/home/livingstonb/GitHub/Continuous_Time_HA/output/two_asset/table.xlsx';

addpath([codedir,'code']);

%% read baseline one-asset for decomposition
fpath = [matdir,'output_oneasset.mat'];
if exist(fpath,'file')
    oneasset = load(fpath);
else
    disp('no one asset baseline found')
    oneasset = [];
end

%% Read .mat files into a cell array

ind = 0;
for run = 1:999
    runstr = num2str(run);
    fpath = [matdir,'output_',runstr,'.mat'];
    if exist(fpath,'file')
        ind = ind+1;
        
        % to save memory
        s{ind} = load(fpath);
        
        if ind > 1
            s{ind}.grd = [];
            s{ind}.KFE = [];
        end
        
        % perform Empc1 - Empc0 decomposition
        decomp_base{ind} = statistics.decomp_baseline(s{1},s{ind});    

        % perform decomp wrt one-asset model
        decomp_oneasset{ind} = statistics.decomp_twoasset_oneasset(oneasset,s{ind});
    else
        continue
    end
end

n = numel(s);

rownames = {'Name'
            'Index'
            'chi0'
            'chi1'
            'chi2'
            'chi1^(-chi2)/(1+chi2)'
            'Rho'
            'Beta (Annualized)'
            'r_a'
            'Illiq Assets'
            'Liq Assets'
            'Total Assets'
            'Wealth Top 10% Share'
            'Wealth Top 1% Share'
            'Gini (Total Assets)'
            's = 0'
            '____HtM STATISTICS'
            'Wealth HtM / Total HtM (at 1/6 qincome)'
            'Wealth HtM / Total HtM (at 1/12 qincome)'
            '____CONSTRAINED IN LIQUID WEALTH'
            'L Wealth == 0'
            'L Wealth <= 0.5% Mean Annual Income'
            'L Wealth <= 1% Mean Annual Income'
            'L Wealth <= 2% * Mean Annual Income'
            'L Wealth <= 5% * Mean Annual Income'
            'L Wealth <= 10% * Mean Annual Income'
            'L Wealth <= 15% * Mean Annual Income'
            'L Wealth <= (1/6) * Own Quarterly Income'
            'L Wealth <= (1/12) * Own Quarterly Income'
            '____CONSTRAINED IN TOTAL WEALTH'
            'Wealth == 0'
            'Wealth <= 0.5% Mean Annual Income'
            'Wealth <= 1% Mean Annual Income'
            'Wealth <= 2% * Mean Annual Income'
            'Wealth <= 5% * Mean Annual Income'
            'Wealth <= 10% * Mean Annual Income'
            'Wealth <= 15% * Mean Annual Income'
            'Wealth <= (1/6) * Own Quarterly Income'
            'Wealth <= (1/12) * Own Quarterly Income'
            '____ILLIQUID ASSETS ADJUSTMENT'
            'Mean |d|/max(a,a_lb)'
            'Median |d|/max(a,a_lb)'
            'Mean chi/|d|, given |d|>0'
            'Median chi/|d|, given |d|>0'
            'Mean chi'
            'Fraction with d == 0'
            '____WEALTH PERCENTILES'
            'Liquid wealth, 10th'
            'Liquid wealth, 25th'
            'Liquid wealth, 50th'
            'Liquid wealth, 90th'
            'Liquid wealth, 99th'
            'Liquid wealth, 99.9th'
            'Total wealth, 10th'
            'Total wealth, 25th'
            'Total wealth, 50th'
            'Total wealth, 90th'
            'Total wealth, 99th'
            'Total wealth, 99.9th'
            '____AVG MPC OUT OF -1e-5 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = -1e-5'
            'QUARTER 2 MPC, shock = -1e-5'
            'QUARTER 3 MPC shock = -1e-5'
            'QUARTER 4 MPC shock = -1e-5'
            'ANNUAL MPC, shock = -1e-5'
            '____AVG MPC OUT OF -0.01 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = -0.01'
            'QUARTER 2 MPC, shock = -0.01'
            'QUARTER 3 MPC shock = -0.01'
            'QUARTER 4 MPC shock = -0.01'
            'ANNUAL MPC, shock = -0.01'
            '____AVG MPC OUT OF -0.1 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = -0.1'
            'QUARTER 2 MPC, shock = -0.1'
            'QUARTER 3 MPC shock = -0.1'
            'QUARTER 4 MPC shock = -0.1'
            'ANNUAL MPC, shock = -0.1'
            '____AVG MPC OUT OF 1e-5 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = 1e-5'
            'QUARTER 2 MPC, shock = 1e-5'
            'QUARTER 3 MPC shock = 1-5'
            'QUARTER 4 MPC shock = 1e-5'
            'ANNUAL MPC, shock = 1e-5'
            '____AVG MPC OUT OF 0.01 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = 0.01'
            'QUARTER 2 MPC, shock = 0.01'
            'QUARTER 3 MPC shock = 0.01'
            'QUARTER 4 MPC shock = 0.01'
            'ANNUAL MPC, shock = 0.01'
            '____AVG MPC OUT OF 0.1 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = 0.1'
            'QUARTER 2 MPC, shock = 0.1'
            'QUARTER 3 MPC shock = 0.1'
            'QUARTER 4 MPC shock = 0.1'
            'ANNUAL MPC, shock = 0.1'
            '____AVG QUARTERLY MPC OUT OF NEWS'
            'QUARTER 1 MPC, shock = -0.01 next quarter'
            'QUARTER 1 MPC, shock = -0.01 next year'
            'ANNUAL MPC, shock = -0.01 next year'
            'QUARTER 1 MPC, shock = 0.01 next quarter'
            'QUARTER 1 MPC, shock = 0.1 next quarter'
            'QUARTER 1 MPC, shock = 0.01 next year'
            'QUARTER 1 MPC, shock = 0.1 next year'
            'ANNUAL MPC, shock = 0.01 next year'
            'ANNUAL MPC, shock = 0.1 next year'
            '___DECOMP OF EM1 AROUND 0'
            'Decomp around 0, RA MPC'
            'Decomp around 0, HtM Effect'
            'Decomp around 0, Non-HtM, constraint'
            'Decomp around 0, Non-HtM, inc risk'
            '___DECOMP OF EM1 AROUND 0.01'
            'Decomp around 0.01, RA MPC'
            'Decomp around 0.01, HtM Effect'
            'Decomp around 0.01, Non-HtM, constraint'
            'decomp Around 0.01, Non-HtM, inc risk'
            '___DECOMP OF EM1 AROUND 0.05'
            'Decomp around 0.05, RA MPC'
            'Decomp around 0.05, HtM Effect'
            'Decomp around 0.05, Non-HtM, constraint'
            'Decomp around 0.05, Non-HtM, inc risk'
            '___DECOMP OF EM1-EM0 (MPC OUT OF 0.01)'
            'Em1 - Em0'
            'Em1 - Em0, Effect of MPC fcn'
            'Em1 - Em0, Effect of distr'
            'Em1 - Em0, Interaction'
            'Distr effect around 0, PHtM only'
            'Distr effect around 0, WHtM only'
            'Distr effect around 0, non-HtM only'
            'Distr effect around 0.01, PHtM only'
            'Distr effect around 0.01, WHtM only'
            'Distr effect around 0.01, non-HtM only'
            'Distr effect around 0.05, PHtM only'
            'Distr effect around 0.05, WHtM only'
            'Distr effect around 0.05, non-HtM only'
            '___DECOMP WRT ONE-ASSET MODEL'
            'E[mpc_twoasset] - E[mpc_oneasset]'
            'Effect of MPC fn'
            'Effect of MPC fn, PHtM term around 0'
            'Effect of MPC fn, WHtM term around 0'
            'Effect of MPC fn, NHtM term around 0'
            'Effect of MPC fn, PHtM term around 0.01'
            'Effect of MPC fn, WHtM term around 0.01'
            'Effect of MPC fn, NHtM term around 0.01'
            'Effect of distr of net worth'
            'Effect of distr of net worth, PHtM term around 0'
            'Effect of distr of net worth, n>0 term around 0'
            'Effect of distr of net worth, PHtM term around 0.01'
            'Effect of distr of net worth, n>0 term around 0.01'
            'Interaction'
            '___DECOMP OF EM1-RA MPC'
            'Em1 - EmRA'
            'Em1 - EmRA, Effect of MPC fcn'
            'Em1 - EmRA, Effect of distr'
            'Em1 - EmRA. Interaction'
            };
        
colnames = {};

T_array = NaN(numel(rownames)-1,n);
for run = 1:n
    colnames{end+1} = s{run}.p.name;
    RHS =  [            run
                        s{run}.p.chi0
                        s{run}.p.chi1
                        s{run}.p.chi2
                        s{run}.stats.adjcosts.chivar
                        s{run}.stats.rho
                        s{run}.stats.beta
                        s{run}.p.r_a
                        s{run}.stats.illiqw
                        s{run}.stats.liqw
                        s{run}.stats.totw
                        s{run}.stats.top10share
                        s{run}.stats.top1share
                        s{run}.stats.wgini
                        s{run}.stats.sav0
                        NaN
                        1 - s{run}.stats.HtM_one_sixth_Q_twealth/s{run}.stats.HtM_one_sixth_Q_lwealth
                        1 - s{run}.stats.HtM_one_twelfth_Q_twealth/s{run}.stats.HtM_one_twelfth_Q_lwealth
                        NaN
                        s{run}.stats.constrained_liq(:)
                        s{run}.stats.HtM_one_sixth_Q_lwealth
                        s{run}.stats.HtM_one_twelfth_Q_lwealth
                        NaN
                        s{run}.stats.constrained(:)
                        s{run}.stats.HtM_one_sixth_Q_twealth
                        s{run}.stats.HtM_one_twelfth_Q_twealth
                        NaN
                        s{run}.stats.adjcosts.mean_d_div_a
                        s{run}.stats.adjcosts.median_d_div_a
                        s{run}.stats.adjcosts.mean_chi_div_d
                        s{run}.stats.adjcosts.median_chi_div_d
                        s{run}.stats.adjcosts.mean_chi
                        s{run}.stats.adjcosts.d0
                        NaN
                        s{run}.stats.lwpercentile(:)
                        s{run}.stats.wpercentile(:)
                        NaN
                        s{run}.stats.mpcs(1).avg_0_quarterly(:)
                        s{run}.stats.mpcs(1).avg_0_annual
                        NaN
                        s{run}.stats.mpcs(2).avg_0_quarterly(:)
                        s{run}.stats.mpcs(2).avg_0_annual
                        NaN
                        s{run}.stats.mpcs(3).avg_0_quarterly(:)
                        s{run}.stats.mpcs(3).avg_0_annual
                        NaN
                        s{run}.stats.mpcs(4).avg_0_quarterly(:)
                        s{run}.stats.mpcs(4).avg_0_annual
                        NaN
                        s{run}.stats.mpcs(5).avg_0_quarterly(:)
                        s{run}.stats.mpcs(5).avg_0_annual
                        NaN
                        s{run}.stats.mpcs(6).avg_0_quarterly(:)
                        s{run}.stats.mpcs(6).avg_0_annual
                        NaN
                        s{run}.stats.mpcs(2).avg_1_quarterly
                        s{run}.stats.mpcs(2).avg_4_quarterly(1)
                        s{run}.stats.mpcs(2).avg_4_annual
                        s{run}.stats.mpcs(5).avg_1_quarterly
                        s{run}.stats.mpcs(6).avg_1_quarterly
                        s{run}.stats.mpcs(5).avg_4_quarterly(1)
                        s{run}.stats.mpcs(6).avg_4_quarterly(1)
                        s{run}.stats.mpcs(5).avg_4_annual
                        s{run}.stats.mpcs(6).avg_4_annual
                        NaN
                        s{run}.stats.decomp_norisk.term1(1)
                        s{run}.stats.decomp_norisk.term2(1)
                        s{run}.stats.decomp_norisk.term3(1)
                        s{run}.stats.decomp_norisk.term4(1)
                        NaN
                        s{run}.stats.decomp_norisk.term1(2)
                        s{run}.stats.decomp_norisk.term2(2)
                        s{run}.stats.decomp_norisk.term3(2)
                        s{run}.stats.decomp_norisk.term4(2)
                        NaN
                        s{run}.stats.decomp_norisk.term1(3)
                        s{run}.stats.decomp_norisk.term2(3)
                        s{run}.stats.decomp_norisk.term3(3)
                        s{run}.stats.decomp_norisk.term4(3)
                        NaN
                        decomp_base{run}.Em1_less_Em0
                        decomp_base{run}.term1
                        decomp_base{run}.term2
                        decomp_base{run}.term3
                        decomp_base{run}.term2a(1)
                        decomp_base{run}.term2b(1)
                        decomp_base{run}.term2c(1)
                        decomp_base{run}.term2a(2)
                        decomp_base{run}.term2b(2)
                        decomp_base{run}.term2c(2)
                        decomp_base{run}.term2a(3)
                        decomp_base{run}.term2b(3)
                        decomp_base{run}.term2c(3)
                        NaN
                        decomp_oneasset{run}.Em1_minus_Em0
                        decomp_oneasset{run}.term1_mpcfn
                        decomp_oneasset{run}.term1a(1)
                        decomp_oneasset{run}.term1b(1)
                        decomp_oneasset{run}.term1c(1)
                        decomp_oneasset{run}.term1a(2)
                        decomp_oneasset{run}.term1b(2)
                        decomp_oneasset{run}.term1c(2)
                        decomp_oneasset{run}.term2_networth
                        decomp_oneasset{run}.term2a(1)
                        decomp_oneasset{run}.term2b(1)
                        decomp_oneasset{run}.term2a(2)
                        decomp_oneasset{run}.term2b(2)
                        decomp_oneasset{run}.term3_interaction
                        NaN
                        s{run}.stats.decompRA.Em1_less_mRA
                        s{run}.stats.decompRA.term1
                        s{run}.stats.decompRA.term2
                        s{run}.stats.decompRA.term3
                        ];
    T_array(:,run) = RHS;
end

% sort by chi2 and chi1
% [T_array,ind] = sortrows(T_array',[3 2]);
% T_array = T_array';
% colnames = {colnames{ind}};

Tnum = T_array;
T_array = num2cell(T_array);
T_array = [colnames;T_array];
T = array2table(T_array);
T.Properties.RowNames = rownames;

writetable(T,xlxpath,'WriteRowNames',true)