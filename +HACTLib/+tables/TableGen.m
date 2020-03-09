classdef TableGen < handle

	properties (SetAccess = protected)
		mpcs_present = false;
		illiquid_mpcs_present = false;
		mpcs_news_present = false;
		include_two_asset_stats = false;
		decomp_norisk_present = false;
		selected_cases;

		n_cols;
	end

	properties
		outdir;
		output = table();
		decomp_baseline;
		decomp_incrisk_alt;
	end

	properties (Abstract)
		fname;
		included_groups;
	end

	methods
		function obj = TableGen(params, stats, use_all)
			obj.outdir = params(1).out_dir;

			if nargin < 3
				use_all = true;
			end
			obj.filter_experiments(params, use_all);

			obj.set_options(params, stats);
			obj.outdir = params(1).out_dir;
		end

		function set_options(obj, params, stats)
			obj.mpcs_present = any([params.ComputeMPCS]);
			obj.illiquid_mpcs_present = any([params.ComputeMPCS_illiquid]);
			obj.mpcs_news_present = any([params.ComputeMPCS_news]);
			obj.include_two_asset_stats = any(~[params.OneAsset]);

			tmp = [stats.decomp_norisk];
			obj.decomp_norisk_present = any([tmp.completed]);

			obj.n_cols = numel(params);
		end

		function filter_experiments(obj, params, use_all)
			if use_all || isempty(obj.included_groups)
				obj.selected_cases = 1:numel(params);
			else
				all_names = {params.group};
				inames = [];
				for ii = 1:numel(all_names)
					if ismember(obj.included_groups, all_names{ii})
						inames = [inames, ii];
					end
				end

				obj.selected_cases = unique(inames);
			end
		end

		function save_table(obj, fname)
			if nargin == 1
				fpath = fullfile(obj.outdir, obj.fname);
			else
				fpath = fullfile(obj.outdir, fname);
			end
			writetable(obj.output, fpath, 'WriteRowNames', true,...
				'WriteVariableNames', false);
		end
	end

	methods (Static)
		function new_table = new_table_with_header(header_name)
			header_formatted = strcat('____', header_name);
			new_table = table({NaN},...
				'VariableNames', {'results'},...
				'RowNames', {header_formatted});
		end

		function output_table = append_to_table(input_table, new_entries)
			table_to_append = table(new_entries(:,2),...
				'VariableNames', {'results'},...
				'RowNames', new_entries(:,1));
			output_table = [input_table; table_to_append];
		end
	end
end