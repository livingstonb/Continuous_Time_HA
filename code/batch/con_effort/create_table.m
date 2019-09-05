clearvars -except stats p

%% This script is used to combine one or more variablesX.mat files. Produces a table.
%% Set FROM_MATFILE = false if running right after model, true if running from .mat file
FROM_MATFILE = true;
codedir = '/home/livingstonb/GitHub/Continuous_Time_HA/';
% codedir = '/Users/Brian-laptop/Documents/GitHub/Continuous_ConEffort/';

if FROM_MATFILE
    % User must set basedir and date, where variablesX.mat files
    % are stored in <basedir>/<date>

    %% Select directory of mat files
    matdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort';
    % matdir = '/Users/Brian-laptop/Documents/midway2_output/8_27_19/';
    xlxpath = '/home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort/table.xlsx';

    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [matdir,'output_',runstr,'.mat'];
        if exist(fpath,'file')
            ind = ind+1;
            
            % to save memory
            s{ind} = load(fpath);
            s{ind}.grd = [];
            s{ind}.grdKFE = [];
            s{ind}.KFE = [];
        else
            continue
        end
    end
    
else
    s = cell(1);
    s{1}.p = p;
    s{1}.stats = stats;
end

n = numel(s);

rownames = {'Name'
            'Definition of h'
            'Index'
            'chi0, coeff on |h|'
            'chi1, coeff on |h|^chi2'
            'chi2'
            'penalty1, coeff on |a|^penalty2 for a<0'
            'penalty2'
            'Rho'
            'Beta (Annualized)'
            'r'
            'Mean Wealth'
            'Wealth Top 10% Share'
            'Wealth Top 1% Share'
            'Gini (Assets)'
            'Fraction of HHs with a < 0'
            '____WEALTH CONSTRAINED'
            'Wealth == 0'
            'Wealth <= 0.5% Mean Annual Income'
            'Wealth <= 1% Mean Annual Income'
            'Wealth <= 2% * Mean Annual Income'
            'Wealth <= 5% * Mean Annual Income'
            'Wealth <= 10% * Mean Annual Income'
            'Wealth <= 15% * Mean Annual Income'
            'Wealth <= (1/6) * Own Quarterly Income'
            'Wealth <= (1/12) * Own Quarterly Income'
            '____WEALTH PERCENTILES'
            '10th'
            '25th'
            '50th'
            '90th'
            '99th'
            '99.9th'
            '____CONSUMPTION ADJUSTMENT'
            'Fraction with h = 0'
            'Fraction with h > 0'
            'Fraction with h < 0'
            'Mean |h|'
            'Mean |h| condl on |h|>0'
            '____CONSUMPTION DISTRIBUTION'
            'c, 10th percentile'
            'c, 95th percentile'
            'c, 99th percentile'
            'c, 99.9th percentile'
            '____AVG SIM MPCs OUT OF -0.01 MEAN ANNUAL INC'
            'QUARTER 1 SIM MPC, shock = -0.01'
            'QUARTER 2 SIM MPC, shock = -0.01'
            'QUARTER 3 SIM MPC, shock = -0.01'
            'QUARTER 4 SIM MPC, shock = -0.01'
            'ANNUAL SIM MPC, shock = -0.01'
            '____AVG SIM MPCs OUT OF -0.1 MEAN ANNUAL INC'
            'QUARTER 1 SIM MPC, shock = -0.1'
            'QUARTER 2 SIM MPC, shock = -0.1'
            'QUARTER 3 SIM MPC, shock = -0.1'
            'QUARTER 4 SIM MPC, shock = -0.1'
            'ANNUAL SIM MPC, shock = 0.1'
            '____AVG SIM MPCs OUT OF 0.01 MEAN ANNUAL INC'
            'QUARTER 1 SIM MPC, shock = 0.01'
            'QUARTER 2 SIM MPC, shock = 0.01'
            'QUARTER 3 SIM MPC, shock = 0.01'
            'QUARTER 4 SIM MPC, shock = 0.01'
            'ANNUAL SIM MPC, shock = 0.01'
            '____AVG SIM MPCs OUT OF 0.1 MEAN ANNUAL INC'
            'QUARTER 1 SIM MPC, shock = 0.1'
            'QUARTER 2 SIM MPC, shock = 0.1'
            'QUARTER 3 SIM MPC, shock = 0.1'
            'QUARTER 4 SIM MPC, shock = 0.1'
            'ANNUAL SIM MPC, shock = 0.1'
            '____AVG SIM MPCs OUT OF 0.01, CONDL ON MPC > 0'
            'QUARTER 1 SIM MPC | MPC > 0, shock = 0.01'
            'QUARTER 2 SIM MPC | MPC > 0, shock = 0.01'
            'QUARTER 3 SIM MPC | MPC > 0, shock = 0.01'
            'QUARTER 4 SIM MPC | MPC > 0, shock = 0.01'
            '____AVG SIM MPCs OUT OF 0.1, CONDL ON MPC > 0'
            'QUARTER 1 SIM MPC | MPC > 0, shock = 0.1'
            'QUARTER 2 SIM MPC | MPC > 0, shock = 0.1'
            'QUARTER 3 SIM MPC | MPC > 0, shock = 0.1'
            'QUARTER 4 SIM MPC | MPC > 0, shock = 0.1'
            '____RESPONSE OUT OF IMMEDIATE SHOCK > 0, SHOCK = 0.01'
            'FRACTION w/ QUARTER 1 MPC > 0, shock=0.01'
            'FRACTION w/ QUARTER 2 MPC > 0, shock=0.01'
            'FRACTION w/ QUARTER 3 MPC > 0, shock=0.01'
            'FRACTION w/ QUARTER 4 MPC > 0, shock=0.01'
            'FRACTION w/ ANNUAL MPC > 0, shock=0.01'
            '____RESPONSE OUT OF IMMEDIATE SHOCK > 0, SHOCK = 0.1'
            'FRACTION w/ QUARTER 1 MPC > 0, shock=0.1'
            'FRACTION w/ QUARTER 2 MPC > 0, shock=0.1'
            'FRACTION w/ QUARTER 3 MPC > 0, shock=0.1'
            'FRACTION w/ QUARTER 4 MPC > 0, shock=0.1'
            'FRACTION w/ ANNUAL MPC > 0, shock=0.1'
            '____AVG SIM MPCs OUT OF NEWS (need to be tested), 0.01 SHOCK'
            'QUARTER 1 MPC, SHOCK NEXT QUARTER, shock = 0.01'
            'QUARTER 1 MPC, SHOCK NEXT YEAR, shock = 0.01'
            'QUARTER 2 MPC, SHOCK NEXT YEAR, shock = 0.01'
            'QUARTER 3 MPC, SHOCK NEXT YEAR, shock = 0.01'
            'QUARTER 4 MPC, SHOCK NEXT YEAR, shock = 0.01'
            '____AVG SIM MPCs OUT OF NEWS (need to be tested), 0.1 SHOCK'
            'QUARTER 1 MPC, SHOCK NEXT QUARTER, shock = 0.1'
            'QUARTER 1 MPC, SHOCK NEXT YEAR, shock = 0.1'
            'QUARTER 2 MPC, SHOCK NEXT YEAR, shock = 0.1'
            'QUARTER 3 MPC, SHOCK NEXT YEAR, shock = 0.1'
            'QUARTER 4 MPC, SHOCK NEXT YEAR, shock = 0.1'
            '____AVG FK MPC OUT OF -1e-5 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = -1e-5'
            'QUARTER 2 MPC, shock = -1e-5'
            'QUARTER 3 MPC shock = -1e-5'
            'QUARTER 4 MPC shock = -1e-5'
            'ANNUAL FK MPC, shock = -1e-5'
            '____AVG FK MPC OUT OF -0.01 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = -0.01'
            'QUARTER 2 MPC, shock = -0.01'
            'QUARTER 3 MPC shock = -0.01'
            'QUARTER 4 MPC shock = -0.01'
            'ANNUAL MPC, shock = -0.01'
            '____AVG FK MPC OUT OF -0.1 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = -0.1'
            'QUARTER 2 MPC, shock = -0.1'
            'QUARTER 3 MPC shock = -0.1'
            'QUARTER 4 MPC shock = -0.1'
            'ANNUAL MPC, shock = -0.1'
            '____AVG FK MPC OUT OF 1e-5 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = 1e-5'
            'QUARTER 2 MPC, shock = 1e-5'
            'QUARTER 3 MPC shock = 1-5'
            'QUARTER 4 MPC shock = 1e-5'
            'ANNUAL MPC, shock = 1e-5'
            '____AVG FK MPC OUT OF 0.01 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = 0.01'
            'QUARTER 2 MPC, shock = 0.01'
            'QUARTER 3 MPC shock = 0.01'
            'QUARTER 4 MPC shock = 0.01'
            'ANNUAL MPC, shock = 0.01'
            '____AVG FK MPC OUT OF 0.1 * MEAN ANNUAL INCOME'
            'QUARTER 1 MPC, shock = 0.1'
            'QUARTER 2 MPC, shock = 0.1'
            'QUARTER 3 MPC shock = 0.1'
            'QUARTER 4 MPC shock = 0.1'
            'ANNUAL MPC, shock = 0.1'
            '____AVG FK QUARTERLY MPC OUT OF NEWS'
            'QUARTERLY MPC, shock = 0.01 next quarter'
            'QUARTERLY MPC, shock = 0.1 next quarter'
            'QUARTERLY MPC, shock = 0.01 next year'
            'QUARTERLY MPC, shock = 0.1 next year'
            };
        
colnames = {};
hdefs = {};

T_array = NaN(numel(rownames)-2,n);
for run = 1:n
    colnames{end+1} = s{run}.p.name;
    hdefs{end+1} = s{run}.p.hdef;
    RHS =  [            run
                        s{run}.p.chi0
                        s{run}.p.chi1
                        s{run}.p.chi2
                        s{run}.p.penalty1
                        s{run}.p.penalty2
                        s{run}.stats.rho
                        s{run}.stats.beta_annualized
                        s{run}.p.r
                        s{run}.stats.wealth
                        s{run}.stats.top10share
                        s{run}.stats.top1share
                        s{run}.stats.wgini
                        s{run}.stats.anegative
                        NaN
                        s{run}.stats.constrained(:)
                        s{run}.stats.HtM_one_sixth_Q_wealth
                        s{run}.stats.HtM_one_twelfth_Q_wealth
                        NaN
                        s{run}.stats.wpercentile(:)
            			NaN
            			s{run}.stats.h0
            			s{run}.stats.hpos
            			s{run}.stats.hneg
            			s{run}.stats.mean_absh_total
            			s{run}.stats.mean_absh_ofchangers
                        NaN
                        s{run}.stats.c10
                        s{run}.stats.c95
                        s{run}.stats.c99
                        s{run}.stats.c999
                        NaN
                        s{run}.stats.sim_mpcs(2).avg_0_quarterly(:)
                        s{run}.stats.sim_mpcs(2).avg_0_annual
                        NaN
                        s{run}.stats.sim_mpcs(3).avg_0_quarterly(:)
                        s{run}.stats.sim_mpcs(3).avg_0_annual
                        NaN
                        s{run}.stats.sim_mpcs(5).avg_0_quarterly(:)
                        s{run}.stats.sim_mpcs(5).avg_0_annual
                        NaN
                        s{run}.stats.sim_mpcs(6).avg_0_quarterly(:)
                        s{run}.stats.sim_mpcs(6).avg_0_annual
                        NaN
                        s{run}.stats.sim_mpcs(5).avg_0_responders(:)
                        NaN
                        s{run}.stats.sim_mpcs(6).avg_0_responders(:)
                        NaN
                        s{run}.stats.sim_mpcs(5).responders_0_quarterly(:)
                        s{run}.stats.sim_mpcs(5).responders_0_annual
                        NaN
                        s{run}.stats.sim_mpcs(6).responders_0_quarterly(:)
                        s{run}.stats.sim_mpcs(6).responders_0_annual
                        NaN
                        s{run}.stats.sim_mpcs(5).avg_1_quarterly
                        s{run}.stats.sim_mpcs(5).avg_4_quarterly(:)
                        s{run}.stats.sim_mpcs(5).avg_4_annual
                        NaN
                        s{run}.stats.sim_mpcs(6).avg_1_quarterly
                        s{run}.stats.sim_mpcs(6).avg_4_quarterly(:)
                        s{run}.stats.sim_mpcs(6).avg_4_annual
                        NaN
                        s{run}.stats.mpcs(1).avg_0_t(1)
                        s{run}.stats.mpcs(1).avg_0_t(2)
                        s{run}.stats.mpcs(1).avg_0_t(3)
                        s{run}.stats.mpcs(1).avg_0_t(4)
                        s{run}.stats.mpcs(1).avg_0_0to4
                        NaN
                        s{run}.stats.mpcs(2).avg_0_t(1)
                        s{run}.stats.mpcs(2).avg_0_t(2)
                        s{run}.stats.mpcs(2).avg_0_t(3)
                        s{run}.stats.mpcs(2).avg_0_t(4)
                        s{run}.stats.mpcs(2).avg_0_0to4
                        NaN
                        s{run}.stats.mpcs(3).avg_0_t(1)
                        s{run}.stats.mpcs(3).avg_0_t(2)
                        s{run}.stats.mpcs(3).avg_0_t(3)
                        s{run}.stats.mpcs(3).avg_0_t(4)
                        s{run}.stats.mpcs(3).avg_0_0to4
                        NaN
                        s{run}.stats.mpcs(4).avg_0_t(1)
                        s{run}.stats.mpcs(4).avg_0_t(2)
                        s{run}.stats.mpcs(4).avg_0_t(3)
                        s{run}.stats.mpcs(4).avg_0_t(4)
                        s{run}.stats.mpcs(4).avg_0_0to4
                        NaN
                        s{run}.stats.mpcs(5).avg_0_t(1)
                        s{run}.stats.mpcs(5).avg_0_t(2)
                        s{run}.stats.mpcs(5).avg_0_t(3)
                        s{run}.stats.mpcs(5).avg_0_t(4)
                        s{run}.stats.mpcs(5).avg_0_0to4
                        NaN
                        s{run}.stats.mpcs(6).avg_0_t(1)
                        s{run}.stats.mpcs(6).avg_0_t(2)
                        s{run}.stats.mpcs(6).avg_0_t(3)
                        s{run}.stats.mpcs(6).avg_0_t(4)
                        s{run}.stats.mpcs(6).avg_0_0to4
                        NaN
                        s{run}.stats.mpcs(5).avg_1_t
                        s{run}.stats.mpcs(6).avg_1_t
                        s{run}.stats.mpcs(5).avg_4_t(1)
                        s{run}.stats.mpcs(6).avg_4_t(1)
                        ];
    T_array(:,run) = RHS;
end

% sort by chi2 and chi1
% [T_array,ind] = sortrows(T_array',[4 3 2]);
% T_array = T_array';
% colnames = {colnames{ind}};
% hdefs = {hdefs{ind}};

Tnum = T_array;
T_array = num2cell(T_array);
T_array = [colnames;hdefs;T_array];
T = array2table(T_array);
T.Properties.RowNames = rownames;

writetable(T,xlxpath,'WriteRowNames',true)