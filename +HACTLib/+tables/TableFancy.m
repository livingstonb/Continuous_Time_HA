classdef TableFancy < handle
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

	methods
		function obj = TableFancy(params, stats)
			obj.outdir = params(1).out_dir;

			obj.set_options(params);
			obj.filter_experiments(params);
		end

		function set_options(obj, params)
			obj.mpcs_present = any([params.ComputeMPCS]);
			obj.illiquid_mpcs_present = any([params.ComputeMPCS_illiquid]);
			obj.mpcs_news_present = any([params.ComputeMPCS_news]);
			obj.include_two_asset_stats = any(~[params.OneAsset]);

			obj.n_cols = numel(params);
		end

		function filter_experiments(obj, params, use_all)
			if nargin < 3
				use_all = true;
			end

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

			obj.n_cols = numel(obj.selected_cases);
		end

		function output_table = create(obj, params, stats)
			output_table = table();
			for ip = 1:obj.n_cols
				p_ip = params(ip);
				stats_ip = stats(ip);

				new_column = obj.intro_stats_table(p_ip, stats_ip);

				tmp = obj.income_stats_table();
				new_column = [new_column; tmp];

				tmp = obj.wealth_stats_table(stats_ip);
				new_column = [new_column; tmp];

				tmp = obj.mpc_size_table(stats_ip);
				new_column = [new_column; tmp];

				tmp = obj.mpc_sign_table(stats_ip);
				new_column = [new_column; tmp];

				column_label = sprintf('Specification%d', ip);
				new_column.Properties.VariableNames = {column_label};
				output_table = [output_table, new_column];
			end
		end

		function out = intro_stats_table(obj, p, stats)
			out = table({p.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_entries = {
				stats.mpcs(5).quarterly
				stats.mpcs(5).annual
				stats.beta_A
			};

			out = obj.construct_from_stats(out, new_entries);
		end

		function out = income_stats_table(obj)
			panel_name = 'Income Statistics';
			out = new_table_with_header(panel_name);

			mean_ann_earnings.value = 1.0;
			mean_ann_earnings.label = 'Mean gross annual earnings';
			mean_ann_earnings.two_asset = false;

			stdev_log_gross_earnings.value = 0.710;
			stdev_log_gross_earnings.label = 'Stdev log ann gross earnings';
			stdev_log_gross_earnings.two_asset = false;

			stdev_log_net_earnings.value = 0.710;
			stdev_log_net_earnings.label = 'Stdev log ann net earnings';
			stdev_log_net_earnings.two_asset = false;

			new_entries = {
				mean_ann_earnings
				stdev_log_gross_earnings
				stdev_log_net_earnings
			};

			out = obj.construct_from_stats(out, new_entries);
		end

		function out = wealth_stats_table(obj, stats)
			panel_name = 'Wealth Statistics';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.median_totw
				stats.median_liqw
				stats.sav0
				stats.constrained_liq{1}
				stats.constrained_liq{2}
				stats.constrained_liq{3}
				stats.constrained_liq{4}
				stats.constrained_liq{5}
				stats.constrained_liq{6}
				stats.liqw_lt_ysixth
				stats.liqw_lt_ytwelfth
				stats.WHtM_over_HtM_biweekly
				stats.WHtM_over_HtM_weekly
				stats.w_top10share
				stats.w_top1share
				stats.wgini
			};

			out = obj.construct_from_stats(out, new_entries);
		end

		function out = mpc_size_table(obj, stats)
			panel_name = 'MPC size effects';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(4).quarterly
				stats.mpcs(6).quarterly
			};

			out = obj.construct_from_stats(out, new_entries);
		end

		function out = mpc_sign_table(obj, stats)
			panel_name = 'MPC sign effects';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(1).quarterly
				stats.mpcs(2).quarterly
				stats.mpcs(3).quarterly
			};

			out = obj.construct_from_stats(out, new_entries);
		end

		function table_out = construct_from_stats(obj, table_in, stats_in)
			vals = {};
			labels = {};
			jj = 1;
			for ii = 1:numel(stats_in)
				if (~obj.include_two_asset_stats) && stats_in{ii}.two_asset
					continue
				end

				vals{jj} = stats_in{ii}.value;
				labels{jj} = stats_in{ii}.label;
				jj = jj + 1;
			end

			table_to_append = table(vals(:),...
				'VariableNames', {'results'},...
				'RowNames', labels(:));

			table_out = [table_in; table_to_append];
		end
	end	
end

function new_table = new_table_with_header(header_name)
	header_formatted = strcat('____', header_name);
	new_table = table({NaN},...
		'VariableNames', {'results'},...
		'RowNames', {header_formatted});
end