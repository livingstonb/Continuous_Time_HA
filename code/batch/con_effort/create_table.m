clearvars -except stats p

%% This script is used to combine one or more variablesX.mat files. Produces a table.
%% Set FROM_MATFILE = false if running right after model, true if running from .mat file
FROM_MATFILE = true;
codedir = '/home/livingstonb/GitHub/Continuous_Time_HA/code';
xlxpath = '/home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort/';
largepath = [xlxpath 'large_table.xlsx'];
smallpath = [xlxpath 'small_table.xlsx'];
% xlxpath = '';
% codedir = '/Users/Brian-laptop/Documents/GitHub/Continuous_ConEffort/';

addpath(codedir)

if FROM_MATFILE
    % User must set basedir and date, where variablesX.mat files
    % are stored in <basedir>/<date>

    %% Select directory of mat files
    matdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output/con_effort/';
    % matdir = '/Users/Brian-laptop/Documents/midway2_output/8_27_19/';

    ind = 0;
    for run = 1:999
        runstr = num2str(run);
        fpath = [matdir,'output_',runstr,'.mat'];
        if exist(fpath,'file')
            ind = ind+1;
            
            % to save memory
            s(ind) = load(fpath);
            s(ind).grd = [];
            s(ind).grdKFE = [];
            s(ind).KFE = [];
        else
            continue
        end
    end
    
else
    s = struct();
    s(1).p = p;
    s(1).stats = stats;
end

n = numel(s);      
nans = num2cell(NaN(1,n));        

%% TABLE WITH ALL STATISTICS
largeTable = 		[{'Name'}, aux.get_all_values(s,'p',1,'name')
			        {'Definition of h'}, aux.get_all_values(s,'p',1,'hdef')
			        {'Index'}, num2cell(1:n)
			        {'chi0, coeff on |h|'}, aux.get_all_values(s,'p',1,'chi0')
			        {'chi1, coeff on |h|^chi2'}, aux.get_all_values(s,'p',1,'chi1')
			        {'chi2'}, aux.get_all_values(s,'p',1,'chi2')
			        {'penalty1, coeff on |a|^penalty2 for a<0'}, aux.get_all_values(s,'p',1,'penalty1')
			        {'penalty2'}, aux.get_all_values(s,'p',1,'penalty2')
			        {'Rho'}, aux.get_all_values(s,'p',1,'rho')
			        {'Beta (Annualized)'}, aux.get_all_values(s,'stats',1,'beta_annualized')
			        {'r'}, aux.get_all_values(s,'p',1,'r_b')
			        {'Mean Wealth'}, aux.get_all_values(s,'stats',1,'wealth')
			        {'Wealth Top 10% Share'}, aux.get_all_values(s,'stats',1,'top10share')
			        {'Wealth Top 1% Share'}, aux.get_all_values(s,'stats',1,'top1share')
			        {'Gini (Assets)'}, aux.get_all_values(s,'stats',1,'wgini')
			        {'____WEALTH CONSTRAINED'}, nans
			        {'Wealth < 0'}, aux.get_all_values(s,'stats',1,'anegative')
			        {'Wealth == 0'}, aux.get_all_values(s,'stats',1,'constrained',1)
			        {'Wealth <= 0.5% Mean Annual Income'}, aux.get_all_values(s,'stats',1,'constrained',2)
			        {'Wealth <= 1% Mean Annual Income'}, aux.get_all_values(s,'stats',1,'constrained',3)
			        {'Wealth <= 2% * Mean Annual Income'}, aux.get_all_values(s,'stats',1,'constrained',4)
			        {'Wealth <= 5% * Mean Annual Income'}, aux.get_all_values(s,'stats',1,'constrained',5)
			        {'Wealth <= 10% * Mean Annual Income'}, aux.get_all_values(s,'stats',1,'constrained',6)
			        {'Wealth <= 15% * Mean Annual Income'}, aux.get_all_values(s,'stats',1,'constrained',7)
			        {'Wealth <= (1/6) * Own Quarterly Income'}, aux.get_all_values(s,'stats',1,'HtM_one_sixth_Q_wealth')
			        {'Wealth <= (1/12) * Own Quarterly Income'}, aux.get_all_values(s,'stats',1,'HtM_one_twelfth_Q_wealth')
			        {'____WEALTH PERCENTILES'}, nans
			        {'10th'}, aux.get_all_values(s,'stats',1,'wpercentile',1)
			        {'25th'}, aux.get_all_values(s,'stats',1,'wpercentile',2)
			        {'50th'}, aux.get_all_values(s,'stats',1,'wpercentile',3)
			        {'90th'}, aux.get_all_values(s,'stats',1,'wpercentile',4)
			        {'99th'}, aux.get_all_values(s,'stats',1,'wpercentile',5)
			        {'99.9th'}, aux.get_all_values(s,'stats',1,'wpercentile',6)
			        {'____CONSUMPTION ADJUSTMENT'}, nans
			        {'Fraction with h = 0'}, aux.get_all_values(s,'stats',1,'h0')
			        {'Fraction with h > 0'}, aux.get_all_values(s,'stats',1,'hpos')
			        {'Fraction with h < 0'}, aux.get_all_values(s,'stats',1,'hneg')
			        {'Mean |h|'}, aux.get_all_values(s,'stats',1,'mean_absh_total')
			        {'Mean |h| condl on |h|>0'}, aux.get_all_values(s,'stats',1,'mean_absh_ofchangers')
			        {'____CONSUMPTION DISTRIBUTION'}, nans
			        {'c, 10th percentile'}, aux.get_all_values(s,'stats',1,'c10')
			        {'c, 95th percentile'}, aux.get_all_values(s,'stats',1,'c95')
			        {'c, 99th percentile'}, aux.get_all_values(s,'stats',1,'c99')
			        {'c, 99.9th percentile'}, aux.get_all_values(s,'stats',1,'c999')
			        {'____AVG SIM MPCs OUT OF -0.01 MEAN ANNUAL INC'}, nans
			        {'QUARTER 1 SIM MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_quarterly',1)
			        {'QUARTER 2 SIM MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_quarterly',2)
			        {'QUARTER 3 SIM MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_quarterly',3)
			        {'QUARTER 4 SIM MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_quarterly',4)
			        {'ANNUAL SIM MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_annual')
			        {'____AVG SIM MPCs OUT OF -0.1 MEAN ANNUAL INC'}, nans
			        {'QUARTER 1 SIM MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_quarterly',1)
			        {'QUARTER 2 SIM MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_quarterly',2)
			        {'QUARTER 3 SIM MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_quarterly',3)
			        {'QUARTER 4 SIM MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_quarterly',4)
			        {'ANNUAL SIM MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_annual')
			        {'____AVG SIM MPCs OUT OF 0.01 MEAN ANNUAL INC'}, nans
			        {'QUARTER 1 SIM MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly',1)
			        {'QUARTER 2 SIM MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly',2)
			        {'QUARTER 3 SIM MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly',3)
			        {'QUARTER 4 SIM MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly',4)
			        {'ANNUAL SIM MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_annual')
			        {'____AVG SIM MPCs OUT OF 0.1 MEAN ANNUAL INC'}, nans
			        {'QUARTER 1 SIM MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',1)
			        {'QUARTER 2 SIM MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',2)
			        {'QUARTER 3 SIM MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',3)
			        {'QUARTER 4 SIM MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',4)
			        {'ANNUAL SIM MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',1)
			        {'____AVG SIM MPCs OUT OF 0.01, CONDL ON MPC > 0'}, nans
			        {'QUARTER 1 SIM MPC | MPC > 0, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly_pos',1)
			        {'QUARTER 2 SIM MPC | MPC > 0, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly_pos',2)
			        {'QUARTER 3 SIM MPC | MPC > 0, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly_pos',3)
			        {'QUARTER 4 SIM MPC | MPC > 0, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly_pos',4)
			        {'____AVG SIM MPCs OUT OF 0.1, CONDL ON MPC > 0'}, nans
			        {'QUARTER 1 SIM MPC | MPC > 0, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly_pos',1)
			        {'QUARTER 2 SIM MPC | MPC > 0, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly_pos',2)
			        {'QUARTER 3 SIM MPC | MPC > 0, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly_pos',3)
			        {'QUARTER 4 SIM MPC | MPC > 0, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly_pos',4)
			        {'____RESPONSE OUT OF IMMEDIATE SHOCK > 0, SHOCK = 0.01'}, nans
			        {'FRACTION w/ QUARTER 1 MPC > 0, shock=0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_0_quarterly',1)
			        {'FRACTION w/ QUARTER 2 MPC > 0, shock=0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_0_quarterly',2)
			        {'FRACTION w/ QUARTER 3 MPC > 0, shock=0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_0_quarterly',3)
			        {'FRACTION w/ QUARTER 4 MPC > 0, shock=0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_0_quarterly',4)
			        {'FRACTION w/ ANNUAL MPC > 0, shock=0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_0_annual')
			        {'____RESPONSE OUT OF IMMEDIATE SHOCK > 0, SHOCK = 0.1'}, nans
			        {'FRACTION w/ QUARTER 1 MPC > 0, shock=0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_0_quarterly',1)
			        {'FRACTION w/ QUARTER 2 MPC > 0, shock=0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_0_quarterly',2)
			        {'FRACTION w/ QUARTER 3 MPC > 0, shock=0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_0_quarterly',3)
			        {'FRACTION w/ QUARTER 4 MPC > 0, shock=0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_0_quarterly',4)
			        {'FRACTION w/ ANNUAL MPC > 0, shock=0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_0_annual')
			        {'____AVG SIM MPCs OUT OF NEWS (need to be tested), 0.01 SHOCK'}, nans
			        {'QUARTER 1 MPC, SHOCK NEXT QUARTER, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_1_quarterly',1)
			        {'QUARTER 1 MPC, SHOCK NEXT YEAR, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_4_quarterly',1)
			        {'QUARTER 2 MPC, SHOCK NEXT YEAR, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_4_quarterly',2)
			        {'QUARTER 3 MPC, SHOCK NEXT YEAR, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_4_quarterly',3)
			        {'QUARTER 4 MPC, SHOCK NEXT YEAR, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_4_quarterly',4)
			        {'ANNUAL MPC, SHOCK NEXT YEAR, shock = 0.01'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_4_annual')
			        {'____AVG SIM MPCs OUT OF NEWS (need to be tested), 0.1 SHOCK'}, nans
			        {'QUARTER 1 MPC, SHOCK NEXT QUARTER, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_1_quarterly',1)
			        {'QUARTER 1 MPC, SHOCK NEXT YEAR, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_4_quarterly',1)
			        {'QUARTER 2 MPC, SHOCK NEXT YEAR, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_4_quarterly',2)
			        {'QUARTER 3 MPC, SHOCK NEXT YEAR, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_4_quarterly',3)
			        {'QUARTER 4 MPC, SHOCK NEXT YEAR, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_4_quarterly',4)
			        {'ANNUAL MPC, SHOCK NEXT YEAR, shock = 0.1'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_4_annual')
			        {'____AVG FK MPC OUT OF -1e-5 * MEAN ANNUAL INCOME'}, nans
			        {'QUARTER 1 MPC, shock = -1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',1)
			        {'QUARTER 2 MPC, shock = -1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',2)
			        {'QUARTER 3 MPC, shock = -1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',3)
			        {'QUARTER 4 MPC, shock = -1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',1,'avg_0_quarterly',4)
			        {'ANNUAL FK MPC, shock = -1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',1,'avg_0_annual')
			        {'____AVG FK MPC OUT OF -0.01 * MEAN ANNUAL INCOME'}, nans
			        {'QUARTER 1 MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',1)
			        {'QUARTER 2 MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',2)
			        {'QUARTER 3 MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',3)
			        {'QUARTER 4 MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'mpcs',2,'avg_0_quarterly',4)
			        {'ANNUAL MPC, shock = -0.01'}, aux.get_all_values(s,'stats',1,'mpcs',2,'avg_0_annual')
			        {'____AVG FK MPC OUT OF -0.1 * MEAN ANNUAL INCOME'}, nans
			        {'QUARTER 1 MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',1)
			        {'QUARTER 2 MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',2)
			        {'QUARTER 3 MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',3)
			        {'QUARTER 4 MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'mpcs',3,'avg_0_quarterly',4)
			        {'ANNUAL MPC, shock = -0.1'}, aux.get_all_values(s,'stats',1,'mpcs',3,'avg_0_annual')
			        {'____AVG FK MPC OUT OF 1e-5 * MEAN ANNUAL INCOME'}, nans
			        {'QUARTER 1 MPC, shock = 1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',1)
			        {'QUARTER 2 MPC, shock = 1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',2)
			        {'QUARTER 3 MPC, shock = 1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',3)
			        {'QUARTER 4 MPC, shock = 1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',4,'avg_0_quarterly',4)
			        {'ANNUAL MPC, shock = 1e-5'}, aux.get_all_values(s,'stats',1,'mpcs',4,'avg_0_annual')
			        {'____AVG FK MPC OUT OF 0.01 * MEAN ANNUAL INCOME'}, nans
			        {'QUARTER 1 MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',1)
			        {'QUARTER 2 MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',2)
			        {'QUARTER 3 MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',3)
			        {'QUARTER 4 MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_0_quarterly',4)
			        {'ANNUAL MPC, shock = 0.01'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_0_annual')
			        {'____AVG FK MPC OUT OF 0.1 * MEAN ANNUAL INCOME'}, nans
			        {'QUARTER 1 MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',1)
			        {'QUARTER 2 MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',2)
			        {'QUARTER 3 MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',3)
			        {'QUARTER 4 MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_0_quarterly',4)
			        {'ANNUAL MPC, shock = 0.1'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_0_annual')
			        {'____AVG FK QUARTERLY MPC OUT OF NEWS'}, nans
			        {'QUARTERLY MPC, shock = 0.01 next quarter'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_1_quarterly',1)
			        {'QUARTERLY MPC, shock = 0.1 next quarter'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_1_quarterly',1)
			        {'QUARTERLY MPC, shock = 0.01 next year'}, aux.get_all_values(s,'stats',1,'mpcs',5,'avg_4_quarterly',1)
			        {'QUARTERLY MPC, shock = 0.1 next year'}, aux.get_all_values(s,'stats',1,'mpcs',6,'avg_4_quarterly',1)];

smallTable = [	{'E[MPC] (-0.01)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_quarterly',1)
				{'E[MPC] (-0.1)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_quarterly',1)
				{'E[MPC] (0.01)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly',1)
				{'E[MPC] (0.1)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly',1)
				{'E[MPC] (0.01 next quarter)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_1_quarterly',1)
				{'E[MPC] (0.1 next quarter)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_1_quarterly',1)
				{'E[MPC|MPC>0] (-0.01)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'avg_0_quarterly_pos',1)
				{'E[MPC|MPC>0] (-0.1)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'avg_0_quarterly_pos',1)
				{'E[MPC|MPC>0] (0.01)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_0_quarterly_pos',1)
				{'E[MPC|MPC>0] (0.1)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_0_quarterly_pos',1)
				{'E[MPC|MPC>0] (0.01 next quarter)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'avg_1_quarterly_pos',1)
				{'E[MPC|MPC>0] (0.1 next quarter)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'avg_1_quarterly_pos',1)
				{'P(MPC>0) (-0.01)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',2,'responders_0_quarterly',1)
				{'P(MPC>0) (-0.1)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',3,'responders_0_quarterly',1)
				{'P(MPC>0) (0.01)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_0_quarterly',1)
				{'P(MPC>0) (0.1)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_0_quarterly',1)
				{'P(MPC>0) (0.01 next quarter)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',5,'responders_1_quarterly',1)
				{'P(MPC>0) (0.1 next quarter)'}, aux.get_all_values(s,'stats',1,'sim_mpcs',6,'responders_1_quarterly',1)
				];

largeTable = largeTable;
smallTable = smallTable';


Tlarge = cell2table(largeTable(:,2:end),'RowNames',largeTable(:,1));

rownames = aux.get_all_values(s,'p',1,'name');
rownames = [{'specification'}, rownames];
Tsmall = cell2table(smallTable,'RowNames',rownames);

if ~isempty(xlxpath) && exist(xlxpath,'dir')
	writetable(Tlarge,largepath,'WriteRowNames',true,'WriteVariableNames',false)
	writetable(Tsmall,smallpath,'WriteRowNames',true,'WriteVariableNames',false)
end
