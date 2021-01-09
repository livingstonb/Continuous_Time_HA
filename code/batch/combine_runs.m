clear

[~, currdir] = fileparts(pwd());
if ~strcmp(currdir, 'Continuous_Time_HA')
    msg = 'The user must cd into the Continuous_Time_HA directory';
    bad_dir = MException('Continuous_Time_HA:master', msg);
    throw(bad_dir);
end

addpath('code');

%% Read .mat files into a cell array
try
    ind = 0;
    for irun = 1:999
        fname = sprintf('output_%d.mat', irun);
        fpath = fullfile('output', 'completed', fname);
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
    output_table = tobj.create(params, stats);

    csvpath = fullfile('output', 'output_table.csv');
    writetable(output_table, csvpath, 'WriteRowNames', true);

    xlxpath = fullfile('output', 'output_table.xlsx');
    writetable(output_table, xlxpath, 'WriteRowNames', true);
catch ME
    HACTLib.aux.display_exception_stack(ME);
    rethrow(ME);
end

if ~isempty(getenv('SLURM_ARRAY_TASK_ID'))
    exit
end