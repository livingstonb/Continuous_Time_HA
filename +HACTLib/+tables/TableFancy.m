classdef TableFancy < handle
	properties (SetAccess = protected)
		mpcs_present = false;
		illiquid_mpcs_present = false;
		mpcs_news_present = false;
		one_asset_only = false;
		two_asset_only = false;
		decomp_baseline_present = false;
		decomp_norisk_present = false;
		selected_cases;

		n_cols;

		current_column;
	end

	properties
		outdir;
		output;
		decomp_incrisk_alt;
	end

	methods
		function obj = TableFancy(params, stats)
			obj.outdir = params(1).out_dir;

			obj.set_options(params, stats);
			obj.filter_experiments(params);
		end

		function set_options(obj, params, stats)
			obj.mpcs_present = any([params.ComputeMPCS]);
			obj.illiquid_mpcs_present = any([params.ComputeMPCS_illiquid]);
			obj.mpcs_news_present = any([params.ComputeMPCS_news]);
			obj.decomp_norisk_present = any(...
				cellfun(@(x) x.decomp_norisk_completed, stats));

			for ii = 1:numel(stats)
				if stats{ii}.decomp_baseline_present
					obj.decomp_baseline_present = true;
					break
				end
			end
			obj.one_asset_only = all([params.OneAsset]);
			obj.two_asset_only = all(~[params.OneAsset]);

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

			obj.output = table();
			for ip = 1:obj.n_cols
				p_ip = params(ip);
				stats_ip = stats{ip};

				obj.current_column = table();
				obj.intro_stats_table(p_ip, stats_ip);

				obj.income_stats_table();
				obj.wealth_stats_table(stats_ip);
				obj.mpc_size_table(stats_ip);
				obj.mpc_sign_table(stats_ip);

				obj.decomp_norisk_table(p_ip, stats_ip);

				shock = stats_ip.mpcs(5).shock.value;
				obj.mpc_comparison(stats_ip, shock);
				obj.mpc_comparison_pct(stats_ip, shock);

				if obj.one_asset_only
					assets = {'b'};
				else
					assets = {'w', 'b', 'a'};
				end
				obj.percentiles_tables(stats_ip, assets);
                obj.illiquid_mpcs_table(stats_ip);
                obj.adj_costs_table(stats_ip);
                obj.other_params_table(stats_ip);
                obj.other_stats_table(stats_ip);

                obj.add_column(ip)
			end

			output_table = obj.output;
		end

		function add_column(obj, ip)
			column_label = sprintf('Specification%d', ip);
			obj.current_column.Properties.VariableNames = {column_label};
			obj.output = [obj.output, obj.current_column];
		end

		function intro_stats_table(obj, p, stats)
			out = table({p.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_entries = {
				stats.mpcs(5).quarterly
				stats.mpcs(5).annual
				stats.beta_A
			};

			obj.update_current_column(out, new_entries);
		end

		function grid_size_table(obj, stats)
			panel_name = 'Grid parameters';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.params.nb
				stats.params.na
				stats.params.b_curv
				stats.params.a_curv
				stats.params.amax
				stats.params.bmax
			};

			obj.update_current_column(out, new_entries);
		end

		function income_stats_table(obj)
			panel_name = 'Income Statistics';
			out = new_table_with_header(panel_name);

			mean_ann_earnings.value = 1.0;
			mean_ann_earnings.label = 'Mean gross annual earnings';
			mean_ann_earnings.indicator = 0;

			stdev_log_gross_earnings.value = 0.710;
			stdev_log_gross_earnings.label = 'Stdev log ann gross earnings';
			stdev_log_gross_earnings.indicator = 0;

			stdev_log_net_earnings.value = 0.710;
			stdev_log_net_earnings.label = 'Stdev log ann net earnings';
			stdev_log_net_earnings.indicator = 0;

			new_entries = {
				mean_ann_earnings
				stdev_log_gross_earnings
				stdev_log_net_earnings
			};

			obj.update_current_column(out, new_entries);
		end

		function wealth_stats_table(obj, stats)
			panel_name = 'Wealth Statistics';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.totw
				stats.liqw
				stats.median_totw
				stats.median_liqw
				stats.sav0
				stats.constrained_liq_pct{1}
				stats.constrained_liq_pct{2}
				stats.constrained_liq_pct{3}
				stats.constrained_liq_pct{4}
				stats.constrained_liq_pct{5}
				stats.constrained_liq_pct{6}
				stats.constrained_liq_dollars{1}
				stats.constrained_liq_dollars{2}
				stats.constrained_liq_dollars{3}
				stats.constrained_liq_dollars{4}
				stats.constrained_liq_dollars{5}
				stats.constrained_liq_dollars{6}
				stats.liqw_lt_ysixth
				stats.liqw_lt_ytwelfth
				stats.WHtM_over_HtM_biweekly
				stats.WHtM_over_HtM_weekly
				stats.w_top10share
				stats.w_top1share
				stats.wgini
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_size_table(obj, stats)
			panel_name = 'MPC size effects';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(4).quarterly
				stats.mpcs(6).quarterly
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_sign_table(obj, stats)
			panel_name = 'MPC sign effects';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.mpcs(1).quarterly
				stats.mpcs(2).quarterly
				stats.mpcs(3).quarterly
			};

			obj.update_current_column(out, new_entries);
		end

		function decomp_norisk_table(obj, p, stats)
			if (obj.two_asset_only) || (~obj.decomp_norisk_present)
				return
			end

			panel_name = 'Decomps of E[MPC] wrt RA and no inc risk, $500 shock';
			out = new_table_with_header(panel_name);

			tmp = stats.mpcs(5).quarterly;
			tmp.label = 'Quarterly MPC (%)'
			new_entries = {
				tmp
				stats.decomp_norisk(1).term1_pct
			};
			obj.update_current_column(out, new_entries);

			for ithresh = 1:numel(p.decomp_thresholds)
				threshold = p.decomp_thresholds(ithresh);
				panel_name = sprintf('For HtM threshold #%d', ithresh);
				out = new_table_with_header(panel_name);

				new_entries = {
					stats.decomp_norisk(ithresh).term2;
					stats.decomp_norisk(ithresh).term3;
					stats.decomp_norisk(ithresh).term4;
				};

				obj.update_current_column(out, new_entries);
			end
		end

		function mpc_comparison(obj, stats, shock)
			if ~obj.decomp_baseline_present
				return
			end

			panel_name = sprintf(...
				'Decomposition of E[MPC] - E[MPC_b] out of %s',...
				shock);
			panel_name = strcat(panel_name,...
				', relative to baseline');
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.decomp_baseline.mean_mpc_diff
				stats.decomp_baseline.term1
				stats.decomp_baseline.term2
				stats.decomp_baseline.term2a(3)
				stats.decomp_baseline.term2b(3)
				stats.decomp_baseline.term2c(3)
				stats.decomp_baseline.term3
			};

			obj.update_current_column(out, new_entries);
		end

		function mpc_comparison_pct(obj, stats, shock)
			if ~obj.decomp_baseline_present
				return
			end

			panel_name = ...
				'Decomposition as % of E[MPC] - E[MPC_b]';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.decomp_baseline.term2_pct
				stats.decomp_baseline.term2a_pct(3)
				stats.decomp_baseline.term2b_pct(3)
				stats.decomp_baseline.term2c_pct(3)
				stats.decomp_baseline.term3_pct
			};

			obj.update_current_column(out, new_entries);
		end

		function percentiles_tables(obj, stats, assets)
			for ia = 1:numel(assets)
				switch assets{ia}
					case 'w'
						panel_name = 'Wealth percentiles';
						new_entries = stats.wpercentiles(:);
					case 'b'
						panel_name = 'Liquid wealth percentiles';
						new_entries = stats.lwpercentiles(:);
					case 'a'
						panel_name = 'Illiquid wealth percentiles';
						new_entries = stats.iwpercentiles(:);
				end
				out = new_table_with_header(panel_name);
				obj.update_current_column(out, new_entries);
			end
		end

		function illiquid_mpcs_table(obj, stats)
			if obj.one_asset_only
				return
			end

			panel_name = 'Illiquid mpcs';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.illiquid_mpcs(1).quarterly
				stats.illiquid_mpcs(2).quarterly
				stats.illiquid_mpcs(3).quarterly
				stats.illiquid_mpcs(4).quarterly
				stats.illiquid_mpcs(5).quarterly
				stats.illiquid_mpcs(6).quarterly
				stats.illiquid_mpcs(4).annual
				stats.illiquid_mpcs(5).annual
				stats.illiquid_mpcs(6).annual
			};

			obj.update_current_column(out, new_entries);
		end

		function adj_costs_table(obj, stats)
			if obj.one_asset_only
				return
			end

			cost_lhs = 'cost(d,a)';
			cost_rhs = 'k0 * |d| + k1 * |d| ^ (1 + k2) / (1 + k2)';
			cost_fn = strcat(cost_lhs, ' = ', cost_rhs);

			panel_name = 'Adjustment cost statistics';
			panel_header = strcat(panel_name, ', ', cost_fn);
			out = new_table_with_header(panel_header);

			new_entries = {
				stats.adjcosts.kappa0
				stats.adjcosts.kappa1
				stats.adjcosts.kappa2
				stats.adjcosts.kappa_var
				stats.adjcosts.a_lb
				stats.adjcosts.mean_cost
				stats.adjcosts.mean_d_div_a
			};

			obj.update_current_column(out, new_entries);
		end

		function other_params_table(obj, stats)
			if obj.one_asset_only
				return
			end

			panel_name = 'Other parameters';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.params.r_a
			};

			obj.update_current_column(out, new_entries);
		end

		function other_stats_table(obj, stats)
			panel_name = 'Other statistics';
			out = new_table_with_header(panel_name);

			new_entries = {
				stats.constrained_illiq_pct{1}
				stats.constrained_illiq_pct{2}
				stats.constrained_illiq_pct{3}
				stats.constrained_illiq_pct{4}
				stats.constrained_illiq_pct{5}
				stats.constrained_illiq_pct{6}
				stats.constrained_illiq_dollars{1}
				stats.constrained_illiq_dollars{2}
				stats.constrained_illiq_dollars{3}
				stats.constrained_illiq_dollars{4}
				stats.constrained_illiq_dollars{5}
				stats.constrained_illiq_dollars{6}
				stats.constrained_pct{1}
				stats.constrained_pct{2}
				stats.constrained_pct{3}
				stats.constrained_pct{4}
				stats.constrained_pct{5}
				stats.constrained_pct{6}
				stats.constrained_dollars{1}
				stats.constrained_dollars{2}
				stats.constrained_dollars{3}
				stats.constrained_dollars{4}
				stats.constrained_dollars{5}
				stats.constrained_dollars{6}
				stats.hhs_paying_wealth_tax
				stats.w_lt_ysixth
				stats.w_lt_ytwelfth
				stats.iwshare_b10
				stats.iwshare_b25
			};

			obj.update_current_column(out, new_entries);
		end

		function update_current_column(obj, table_in, stats_in)
			if nargin < 3
				tmp = table_in;
			else
				tmp = obj.construct_from_stats(table_in, stats_in);
			end
			obj.current_column = [obj.current_column; tmp];
		end

		function table_out = construct_from_stats(obj, table_in, stats_in)
			vals = {};
			labels = {};
			jj = 1;
			for ii = 1:numel(stats_in)
				if obj.one_asset_only && (stats_in{ii}.indicator == 2)
					continue
				elseif obj.two_asset_only && (stats_in{ii}.indicator == 1)
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