clear
load('/home/brian/Documents/temp/stats_vars.mat');

stats = HACTLib.model_objects.Statistics(p, income, grdKFE, KFE);
stats.compute_statistics();