classdef OtherPanels
	methods (Static)
		function out = intro_panel(stats, p, group)
			if isempty(p.label)
				param_label = p.name;
			else
				param_label = p.label;
			end

			if nargin < 3
				group = '';
			end

			out = table({param_label},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			% if strcmp(group, 'Q2')
				
			% end

			new_entries = {	'kappa1, adj cost coeff', true
							'kappa2, adj cost power', true
							'Quarterly MPC (%)', false
				            'Annual MPC (%)', false
				            'Beta (Annualized)', false
				};
			new_entries(:,3) = ...
				{	round(stats.adjcosts.kappa1, 3)
					round(stats.adjcosts.kappa2, 3)
					round(stats.mpcs(5).avg_0_quarterly * 100, 1)
                    round(stats.mpcs(5).avg_0_annual * 100, 1)
                    round(stats.beta_annualized, 3)
				};

			if ~obj.include_two_asset_stats
				new_entries([new_entries{:,2}],:) = [];
			end
			out = tables.TableGen.append_to_table(out,...
				new_entries(:,[1 3]));
		end

		function out = wealth_panel(values, panel_name)
			if nargin == 1
				panel_name = 'Panel B: Wealth statistics';
			end

			out = tables.TableGen.new_table_with_header(panel_name);

			% Mean assets and saving
			new_labels = {	'Mean wealth', false
							'Median wealth', false
							'Mean liquid assets', true
							'Median liquid assets', true
				};
			new_entries(:,3) = ...
				{			stats.totw
							stats.median_totw
							stats.liqw
							stats.median_liqw
				};
			new_entries(:,3) = aux.cellround(new_entries(:,3), 3);
			if ~obj.include_two_asset_stats
				new_entries([new_entries{:,2}],:) = [];
			end
			out = tables.TableGen.append_to_table(out,...
				new_entries(:,[1 3]));

			% Fraction with saving = 0
			new_entries = {'s == 0', round(stats.sav0, 3)};
			out = tables.TableGen.append_to_table(out,...
				new_entries);

			% Fraction with assets or cash<= some value
			new_entries = {	'liq. w. <= 0'
				            'liq. w. <= 0.5% mean ann inc'
				            'liq. w. <= 1% mean ann inc'
				            'liq. w. <= 2% mean ann inc'
				            'liq. w. <= 5% mean ann inc'
				            'liq. w. <= 10% mean ann inc'
				            'liq. w. <= 15% mean ann inc'
				            'liq. w. <= 1/6 own quarterly income'
				            'liq. w. <= 1/12 own quarterly income'
				};
			new_entries(:,2) = ...
				[	num2cell(stats.constrained_liq(:))
					stats.HtM_one_sixth_Q_lwealth
					stats.HtM_one_twelfth_Q_lwealth
				];
			new_entries(:,2) = aux.cellround(new_entries(:,2), 3);
			out = tables.TableGen.append_to_table(out,...
				new_entries);

			% Fraction WHtM / HtM
			if obj.include_two_asset_stats
				new_entries = {	'WHtM/HtM at 1/6 inc', stats.ratio_WHtM_HtM_sixth
					            'WHtM/HtM at 1/12 inc', stats.ratio_WHtM_HtM_twelfth
					};
				new_entries(:,2) = aux.cellround(new_entries(:,2), 3);
				out = tables.TableGen.append_to_table(out,...
					new_entries);
			end

			% Gini
			new_entries = { 'Gini coefficient', round(stats.wgini, 3) };
			out = tables.TableGen.append_to_table(out,...
					new_entries);
		end
	end
end
