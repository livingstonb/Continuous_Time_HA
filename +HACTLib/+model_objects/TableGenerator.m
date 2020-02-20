classdef TableGenerator < handle
	properties (SetAccess = private)
		mpcs_present = false;
		mpcs_news_present = false;
		include_two_asset_stats = false;

		n_cols;
	end

	properties
		decomp_incrisk;
		decomp_repagent;
		decomp_incrisk_alt;
	end

	methods
		function output_table = create(obj, cases)
			
			obj.n_cols = numel(cases);
			[params, stats] = create_structures(cases);

			obj.set_options(params);

			output_table = table();
			for ip = 1:numel(cases)
				p = params(ip);
				sts = stats(ip);
				new_column = obj.intro_stats_table(p, sts);

				if obj.include_two_asset_stats
					tmp = htm_stats_table(p, stats);
					new_column = [new_column; tmp];
				end

				tmp = wconstraint_stats_table(p, stats, 'liquid');
				new_column = [new_column; tmp];

				if obj.include_two_asset_stats
					tmp = wconstraint_stats_table(p, stats, 'total');
					new_column = [new_column; tmp];

					tmp = illiq_adj_table(p, stats);
					new_column = [new_column; tmp];
				end

				tmp = obj.wealth_pct_table(stats);
				new_column = [new_column; tmp];

				output_table = [output_table, new_column];
			end
		end

		function set_options(obj, params)
			obj.mpcs_present = any([params.ComputeMPCS]);
			obj.mpcs_news_present = any([params.ComputeMPCS_news]);
			obj.include_two_asset_stats = any(~[params.OneAsset]);
		end

		function out = intro_stats_table(obj, p, stats)
			out = table({p.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_entries = {	'chi1', p.chi1, true
							'chi2', p.chi2, true
							'chi1^(-chi2)/(1+chi2)', stats.adjcosts.chivar, true
							'rho', p.rho, false
							'beta (annualized)', stats.beta_annualized, false
							'r_a', p.r_a, true
							'Mean Illiq Assets', stats.illiqw, true
							'Mean Liq Assets', stats.liqw, false
							'Mean Total Assets', stats.totw, true
							'Wealth, Top 10% Share', stats.top10share, false
							'Wealth, Top 1% Share', stats.top1share, false
							'Gini (total assets)', stats.wgini, false
							's = 0', stats.sav0, false
				};

			if ~obj.include_two_asset_stats
				new_entries([new_entries{:,3}],:) = [];
			end

			out = append_to_table(out, new_entries);
		end

		function out = wealth_pct_table(obj, stats)
			header_name = 'WEALTH PERCENTILES';
			out = new_table_with_header(header_name);

			new_entries = {	'Liquid wealth, 10th', stats.lwpercentile(1), false
							'Liquid wealth, 25th', stats.lwpercentile(2), false
							'Liquid wealth, 50th', stats.lwpercentile(3), false
							'Liquid wealth, 90th', stats.lwpercentile(4), false
							'Liquid wealth, 99th', stats.lwpercentile(5), false
							'Liquid wealth, 99.9th', stats.lwpercentile(6), false
							'Total wealth, 10th', stats.wpercentile(1), true
							'Total wealth, 25th', stats.wpercentile(2), true
							'Total wealth, 50th', stats.wpercentile(3), true
							'Total wealth, 90th', stats.wpercentile(4), true
							'Total wealth, 99th', stats.wpercentile(5), true
							'Total wealth, 99.9th', stats.wpercentile(6), false
				};

			if ~obj.include_two_asset_stats
				new_entries([new_entries{:,3}],:) = [];
			end

			out = append_to_table(out, new_entries);
		end
	end
end

function [params, stats] = create_structures(cases)
	for ip = 1:numel(cases)
		this_case = cases(ip);

		params(ip) = this_case.p;
		stats(ip) = this_case.stats;
	end
end

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

function output_table = htm_stats_table(p, stats)
	header_name = 'HtM STATISTICS';
	output_table = new_table_with_header(header_name);

	new_entries = {	'Wealthy HtM / Total HtM (at 1/6 qincome)', stats.ratio_WHtM_HtM_sixth
					'Wealthy HtM / Total HtM (at 1/12 qincome)', stats.ratio_WHtM_HtM_twelfth
				};

	output_table = append_to_table(output_table, new_entries);
end

function output_table = wconstraint_stats_table(p, stats, asset)
	if strcmp(asset, 'liquid')
		header_label = 'LIQUID';
		short_label = "L ";
		avar = stats.constrained_liq;
	elseif strcmp(asset, 'total')
		header_label = 'TOTAL';
		short_label = '';
		avar = stats.constrained;
	end
	header_name = sprintf('CONSTRAINED IN %s WEALTH', header_label);
	output_table = new_table_with_header(header_name);

	nconstraints = numel(avar);
	new_entries = cell(nconstraints+2, 2);

	for ii = 1:nconstraints
		constraint = p.epsilon_HtM(ii) * 100;
		if ii == 1
			label = sprintf('%sWealth == 0', short_label);
		else
			label = sprintf('%sWealth <= %g%% Mean Ann Inc', short_label, constraint);
		end
		value = avar(ii);

		new_entries{ii,1} = label;
		new_entries{ii,2} = value;
	end

	new_entries{nconstraints+1,1} = sprintf('%sWealth <= (1/6) Own Quart Inc', short_label);
	new_entries{nconstraints+2,1} = sprintf('%sWealth <= (1/12) Own Quart Inc', short_label);

	if strcmp(asset, 'liquid')
		new_entries{nconstraints+1,2} = stats.HtM_one_sixth_Q_lwealth;
		new_entries{nconstraints+2,2} = stats.HtM_one_twelfth_Q_lwealth;
	else
		new_entries{nconstraints+1,2} = stats.HtM_one_sixth_Q_twealth;
		new_entries{nconstraints+2,2} = stats.HtM_one_twelfth_Q_twealth;
	end

	output_table = append_to_table(output_table, new_entries);
end

function output_table = illiq_adj_table(p, stats)
	header_name = 'ILLIQUID ASSETS ADJUSTMENT';
	output_table = new_table_with_header(header_name);

	new_entries = {	'Mean |d|/max(a,a_lb)', stats.adjcosts.mean_d_div_a
					'Median |d|/max(a,a_lb)', stats.adjcosts.median_d_div_a
					'Mean chi/|d|, given |d|>0', stats.adjcosts.mean_chi_div_d
					'Median chi/|d|, given |d|>0', stats.adjcosts.median_chi_div_d
					'Mean chi', stats.adjcosts.mean_chi
					'Fraction with d == 0', stats.adjcosts.d0
				};

	output_table = append_to_table(output_table, new_entries);
end


