classdef TableGenerator
	properties (SetAccess = private)
		output_table;
		current_column;

		n_cols;

		intro_stats_labels = {
			'Rho'
            'Beta (Annualized)'
            
			};

		intro_stats_labels_illiq = {
			'r_a'
			};
	end

	properties
		decomp_incrisk;
		decomp_repagent;
		decomp_incrisk_alt;
	end

	methods
		function obj = TableGenerator(cases)
			obj.n_cols = numel(cases);

			for ip = 1:n_cols
				this_case = cases(ip);

				params(ip) = this_case.p;
				stats(ip) = this_case.stats;

				if all([params.OneAsset])


			obj.non_header_rows = [
					intro_stats_labels,...
				];
		end

		function output_table = create(obj, cases)
			output_table = table();
			n_cols = numel(cases);

			for ip = 1:n_cols
				this_case = cases(ip);
				params(ip) = this_case.p;
				stats(ip) = this_case.stats;
			end

			if all([params.OneAsset])
				obj.rows_to_drop

			
		end

		function out = intro_stats_table(obj, p, stats)
			out = table({p.name},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_labels = obj.intro_stats_labels;

			new_entries = {  
				};

			obj.add_non_header_rows(new_labels);
			out = append_to_table(out, new_entries, intro_stats_labels);
		end
	end
end

