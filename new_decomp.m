
clear
baseline = load('output/output_301.mat');
no_bc = load('output/HA_no_borr_constraint.mat');

addpath('code');
decomp = statistics.decomp_nobc_norisk(baseline, no_bc);