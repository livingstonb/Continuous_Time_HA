clear

server = false;

if ~server
    basedir = '/home/brian/Documents/GitHub/Continuous_Time_HA';
    matdir = '/home/brian/Documents/GitHub/Continuous_Time_HA/output/server2';
    xlxdir = '/home/brian/Documents/GitHub/Continuous_Time_HA/output';
else
    basedir = '/home/livingstonb/GitHub/Continuous_Time_HA';
    matdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output';
    xlxdir = '/home/livingstonb/GitHub/Continuous_Time_HA/output';
end

addpath(basedir);
addpath([basedir '/code']);

%% Read .mat files into a cell array
ind = 0;
for irun = 1:999
    fname = sprintf('output_%d.mat', irun);
    fpath = fullfile(matdir, fname);
    if exist(fpath,'file')
        ind = ind + 1;
        
        s(ind) = load(fpath);
        if ind > 1
            s(ind).grd = [];
            s(ind).KFE = [];
        end
        
        params(ind) = s(ind).p;
        stats{ind} = s(ind).stats;

        % perform Empc1 - Empc0 decomposition
        decomp_base(ind) = statistics.decomp_baseline(s(1), s(ind));

        % perform decomp wrt one-asset model
        % decomp_oneasset(ind) = statistics.decomp_twoasset_oneasset(oneasset,s(ind));
        
        stats{ind} = HACTLib.aux.add_comparison_decomps(params(ind),...
            stats{ind}, decomp_base(ind));
    end
end

fprintf('%d experiments were found..\n\n', ind)

tobj = HACTLib.tables.StatsTable(params, stats);
output_table = tobj.create(params, stats)


csvpath = fullfile(xlxdir, 'output_table.csv');
writetable(output_table, csvpath, 'WriteRowNames', true);

xlxpath = fullfile(xlxdir, 'output_table.xlsx');
writetable(output_table, xlxpath, 'WriteRowNames', true);