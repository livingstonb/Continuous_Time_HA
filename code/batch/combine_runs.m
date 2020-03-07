clear

server = true;

if ~server
    basedir = '/home/brian/Documents/GitHub/Continuous_Time_HA';
    matdir = '/home/brian/Documents/GitHub/Continuous_Time_HA/output';
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

        % % perform Empc1 - Empc0 decomposition
        % if isequal(s(ind).grdKFE.b.vec, s(1).grdKFE.b.vec)
        %     decomp_base(ind) = statistics.decomp_baseline(s(1),s(ind));
        % else
        %     decomp_base(ind) = statistics.decomp_baseline(s(1),s(1));
        % end

        % perform decomp wrt one-asset model
        % decomp_oneasset(ind) = statistics.decomp_twoasset_oneasset(oneasset,s(ind));
    end
end

table_gen = HACTLib.model_objects.TableGenerator();
output_table = table_gen.create(s);

xlxpath = fullfile(xlxdir, 'output_table.xlsx');
writetable(output_table, xlxpath, 'WriteRowNames', true);