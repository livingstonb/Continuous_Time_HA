classdef Statistics < handle

	properties
		pmf;
		rho;
		beta_Q;
		beta_A;
		illiqw;
		liqw;
		totw;
		sav0;

		mpcs_over_ss;
		pmf_b;
		cdf_b;
		pmf_a;
		cdf_a;
		pmf_b_a;

		lwpercentiles;
		iwpercentiles;
		wpercentiles;

		median_liqw;
		median_illiqw;
		median_totw;

		w_top10share;
		w_top1share;
		lw_top10share;
		lw_top1share;
		iw_top10share
		iw_top1share
		iwshare_b10;
		iwshare_b25;
		wgini;

		constrained;
		constrained_pct;
		constrained_dollars;
		constrained_liq;
		constrained_liq_pct;
		constrained_liq_dollars;
		constrained_illiq;
		constrained_illiq_pct;
		constrained_illiq_dollars;
		hhs_paying_wealth_tax;
		
		w_lt_ysixth;
		w_lt_ytwelfth;
		liqw_lt_ysixth;
		liqw_lt_ytwelfth;

		WHtM_over_HtM_biweekly;
		WHtM_over_HtM_weekly;

		adjcosts;

		mpcs;
		illiquid_mpcs;
		mpcs_news_one_quarter;
		mpcs_news_one_year;

		decomp_norisk_completed;
		decomp_norisk;
		decomp_RA;
		decomp_baseline_present = false;

		params = struct();

		other = struct();
	end

	properties (Access=protected)
		p;
		income;
		grdKFE;
		model;

		wealth_sorted;
		wealthmat;

		pmf_w;
		cdf_w;

		nb;
		na;
        
        kernel_options;

        lw_interp;
        iw_interp;
        w_interp;
	end

	methods
		function obj = Statistics(p, income, grdKFE, model)
			obj.p = p;
			obj.income = income;
			obj.grdKFE = grdKFE;
			obj.model = model;

			obj.na = p.na_KFE;
			obj.nb = p.nb_KFE;

			tmp = grdKFE.b.vec + grdKFE.a.wide;
			obj.wealthmat = tmp;
			obj.wealth_sorted = sortrows(tmp(:));

			obj.pmf = model.g .* grdKFE.trapezoidal.matrix;
		end

		function compute_statistics(obj, kernel_options)
			if nargin == 1
				obj.kernel_options = struct();
			else
				obj.kernel_options = kernel_options;
            end

			obj.add_params();
			obj.compute_intro_stats();
			obj.construct_distributions();
			obj.compute_percentiles();
			obj.compute_inequality();
			obj.compute_constrained();
			obj.compute_deposit_stats();
		end

		function add_params(obj)
			obj.params.nb = sfill(obj.p.nb,...
				'nb, no. of points on liquid grid');
			obj.params.na = sfill(obj.p.na,...
				'na, no. of points on illiquid grid', 2);
			obj.params.b_curv = sfill(obj.p.b_gcurv_pos,...
				'Liquid grid curvature', 2);
			obj.params.a_curv = sfill(obj.p.a_gcurv,...
				'Illiquid grid curvature', 2);
			obj.params.bmax = sfill(obj.p.bmax,...
				'Liquid grid max');
			obj.params.amax = sfill(obj.p.amax,...
				'Illiquid grid max');
			obj.params.r_a = sfill(obj.p.r_a,...
				'Illiquid return (quarterly)');
		end

		function add_mpcs(obj, mpc_obj, simulated)
			if mpc_obj.options.liquid_mpc
				asset_indicator = 0;
				mpc_type = '';
			else
				asset_indicator = 2;
				mpc_type = 'illiq';
			end
			sfill2 = @(x,y) sfill(x, y, asset_indicator);

			empty_stat = sfill2([], []);

			empty_mpc_struct = struct(...
				'shock_normalized', empty_stat,...
				'shock', empty_stat,...
				'quarterly', empty_stat,...
				'annual', empty_stat...
			);

			nshocks = numel(obj.p.mpc_shocks);
			for ishock = 1:nshocks
				shock = obj.p.mpc_shocks(ishock);
				shock_label = obj.p.quantity2label(shock);
				mpcs_stats(ishock) = empty_mpc_struct;

				mpcs_stats(ishock).shock_normalized = sfill2(...
					shock * 100, 'Shock size, (% of mean ann inc)');

				mpcs_stats(ishock).shock = sfill2(shock_label,...
					'Shock size');

				tmp = 100 * mpc_obj.mpcs(ishock).quarterly(1);
				label = sprintf(...
					'Quarterly %s MPC (%%), out of %s',...
					mpc_type, shock_label);
				mpcs_stats(ishock).quarterly = sfill2(tmp, label);

				tmp = 100 * mpc_obj.mpcs(ishock).annual;
				label = sprintf(...
					'Annual %s MPC (%%), out of %s',...
					mpc_type, shock_label);
				mpcs_stats(ishock).annual = sfill2(tmp, label);
			end

			if mpc_obj.options.liquid_mpc
				obj.mpcs = mpcs_stats;

				obj.mpcs_over_ss = cell(1, nshocks);
				for ishock = 1:nshocks
					obj.mpcs_over_ss{ishock} =  mpc_obj.mpcs(ishock).mpcs(:,1);
				end
			else
				obj.illiquid_mpcs = mpcs_stats;
			end
		end

		function add_mpcs_news(obj, mpc_obj)
			empty_stat = sfill([], []);

			% Shock in one quarter
			empty_mpc_struct = struct(...
				'shock_normalized', empty_stat,...
				'shock', empty_stat,...
				'quarterly', empty_stat...
			);

			mpcs_stats = empty_mpc_struct;
			nshocks = numel(obj.p.mpc_shocks);
			for ishock = 1:nshocks
				shock = obj.p.mpc_shocks(ishock);
				shock_label = obj.p.quantity2label(shock);

				mpcs_stats(ishock) = empty_mpc_struct;

				mpcs_stats(ishock).shock_normalized = sfill(...
					shock * 100, 'Size of shock next quarter, (% of mean ann inc)');

				mpcs_stats(ishock).shock = sfill(shock_label,...
					'Size of shock next quarter');

				tmp = 100 * mpc_obj.mpcs(ishock).avg_1_quarterly;
				label = sprintf(...
					'Quarterly MPC (%%), out of %s next quarter',...
					shock_label);
				mpcs_stats(ishock).quarterly = sfill(tmp, label);
			end
			obj.mpcs_news_one_quarter = mpcs_stats;

			% Shock in one year
			empty_mpc_struct = struct(...
				'shock_normalized', empty_stat,...
				'shock', empty_stat,...
				'quarterly', empty_stat,...
				'annual', empty_stat...
			);

			mpcs_stats = empty_mpc_struct;
			nshocks = numel(obj.p.mpc_shocks);
			for ishock = 1:nshocks
				shock = obj.p.mpc_shocks(ishock);
				shock_label = obj.p.quantity2label(shock);

				mpcs_stats(ishock) = empty_mpc_struct;

				mpcs_stats(ishock).shock_normalized = sfill(...
					shock * 100, 'Size of shock next year, (% of mean ann inc)');

				mpcs_stats(ishock).shock = sfill(shock_label,...
					'Size of shock next year');

				tmp = 100 * mpc_obj.mpcs(ishock).avg_4_quarterly;
				label = sprintf(...
					'Quarterly MPC (%%), out of %s next year',...
					shock_label);
				mpcs_stats(ishock).quarterly = sfill(tmp, label);

				tmp = 100 * mpc_obj.mpcs(ishock).avg_4_annual;
				label = sprintf(...
					'Annual MPC (%%), out of %s next year',...
					shock_label);
				mpcs_stats(ishock).annual = sfill(tmp, label);
			end
			obj.mpcs_news_one_year = mpcs_stats;
		end

		function add_decomps(obj, decomp_obj)
			empty_stat = sfill([], []);

			% Decomposition wrt no inc risk model
			empty_decomp_struct = struct(...
				'htm_level', empty_stat,...
				'term1', empty_stat, 'term1_pct', empty_stat,...
				'term2', empty_stat, 'term3', empty_stat,...
				'term4', empty_stat...
			);
			n = numel(obj.p.decomp_thresholds);

			norisk = decomp_obj.results_norisk;
			obj.decomp_norisk = empty_decomp_struct;
			obj.decomp_norisk_completed = norisk.completed;
			for it = 1:n
				lvl = obj.p.decomp_thresholds(it);

				obj.decomp_norisk(it) = empty_decomp_struct;
				obj.decomp_norisk(it).htm_level = sfill(...
					lvl, 'HtM threshold');

				obj.decomp_norisk(it).term1 = sfill(...
					norisk.term1(it), 'RA MPC');

				obj.decomp_norisk(it).term1_pct = sfill(...
					norisk.term1(it) * 100, 'RA MPC (%)');

				tmp = sprintf('HtM effect (abar <= %g)', lvl);
				obj.decomp_norisk(it).term2 = sfill(...
					norisk.term2(it), tmp);

				tmp = sprintf(...
					'Non-HtM (abar > %g), constraint effect ', lvl);
				obj.decomp_norisk(it).term3 = sfill(...
					norisk.term3(it), tmp);

				tmp = sprintf(...
					'Non-HtM (abar > %g), inc risk effect ', lvl);
				obj.decomp_norisk(it).term4 = sfill(...
					norisk.term4(it), tmp);
			end

			% Decomposition wrt RA model
			n = numel(obj.p.decomp_thresholds);

			results_RA = decomp_obj.results_RA;
			obj.decomp_RA = struct();
			obj.decomp_RA.mpc_RA = sfill(decomp_obj.m_ra,...
				'MPC of RA model');
			obj.decomp_RA.mean_mpc_diff = sfill(...
				results_RA.Em1_less_mRA,...
				'E[MPC] - MPC_ra');
			obj.decomp_RA.term1 = sfill(...
				results_RA.term1,...
				'Effect of MPC fcn');
			obj.decomp_RA.term2 = sfill(...
				results_RA.term2,...
				'Effect of distribution');
			obj.decomp_RA.term3 = sfill(...
				results_RA.term3,...
				'Interaction');
		end

		function clean(obj)
			obj.p = [];
			obj.income = [];
			obj.grdKFE = [];
			obj.model = [];
			obj.wealthmat = [];
			obj.wealth_sorted = [];
			obj.pmf_w = [];
			obj.cdf_w = [];
			obj.lw_interp = [];
			obj.iw_interp = [];
			obj.w_interp = [];
		end
	end

	methods (Access=protected)
		function compute_intro_stats(obj)
			obj.rho = sfill(obj.p.rho, 'rho');
		    obj.beta_Q = sfill(exp(-obj.p.rho), 'beta (quarterly)');
		    obj.beta_A = sfill(exp(-4 * obj.p.rho), 'beta (annualized)');

		    tmp = obj.expectation(obj.grdKFE.a.wide);
		    obj.illiqw = sfill(tmp, 'Mean illiquid wealth', 2);

		    tmp = obj.expectation(obj.grdKFE.b.vec);
		    obj.liqw = sfill(tmp, 'Mean liquid wealth');

		    tmp = obj.expectation(...
		    	obj.grdKFE.b.vec + obj.grdKFE.a.wide);
		    obj.totw = sfill(tmp, 'Mean total wealth', 2);

		    tmp = obj.expectation(obj.model.s==0);
		    obj.sav0 = sfill(tmp, 's = 0');
		end

		function construct_distributions(obj)
			import HACTLib.aux.multi_sum

		    [obj.pmf_b, obj.cdf_b] = obj.marginal_dists(1);
		    [obj.pmf_a, obj.cdf_a] = obj.marginal_dists(2);
		    [obj.pmf_w, obj.cdf_w] = obj.marginal_dists([1, 2]);
		    obj.pmf_b_a = multi_sum(obj.pmf, [3, 4]);

		    obj.lw_interp = get_interpolant(obj.grdKFE.b.vec,...
				obj.pmf, [1], obj.kernel_options);
		    obj.iw_interp = get_interpolant(obj.grdKFE.a.vec,...
				obj.pmf, [2], obj.kernel_options);

		    obj.w_interp = get_interpolant(obj.wealthmat,...
				obj.pmf, [1, 2], obj.kernel_options);
		end

		function compute_percentiles(obj)
			npct = numel(obj.p.wpercentiles);
		    obj.lwpercentiles = cell(1, npct);
		    obj.iwpercentiles = cell(1, npct);
		    obj.wpercentiles = cell(1, npct);
			for ip = 1:npct
				pct_at = obj.p.wpercentiles(ip);

				tmp_b = sprintf('b, %gth pctile', pct_at);
				obj.lwpercentiles{ip} = sfill(...
					obj.lw_interp.icdf(pct_at/100), tmp_b);

				tmp_a = sprintf('a, %gth pctile', pct_at);
				obj.iwpercentiles{ip} = sfill(...
					obj.iw_interp.icdf(pct_at/100), tmp_a, 2);

				tmp_w = sprintf('w, %gth pctile', pct_at);
				obj.wpercentiles{ip} = sfill(...
					obj.w_interp.icdf(pct_at/100), tmp_w, 2);
			end

			obj.median_liqw = sfill(obj.lw_interp.icdf(0.5),...
				'b, median');
			obj.median_illiqw = sfill(obj.iw_interp.icdf(0.5),...
				'a, median', 2);
			obj.median_totw = sfill(obj.w_interp.icdf(0.5),...
				'w, median', 2);
		end

		function compute_inequality(obj)
			import HACTLib.aux.interp_integral_alt
			import HACTLib.aux.direct_gini
			import HACTLib.aux.multi_sum

			wshares_kernel_options = obj.kernel_options;
			wshares_kernel_options.x_transform = @(x) log(0.001 + 100*x);
		    wshares_kernel_options.x_revert = @(z) (exp(z) - 0.001) / 100;

			% Top wealth shares
			pmf_w = multi_sum(obj.pmf, [3, 4]);
			tmp = sortrows([obj.wealthmat(:), pmf_w(:)]);
			values_w = cumsum(tmp(:,1) .* tmp(:,2) / obj.totw.value);
			wcumshare_interp = get_interpolant(values_w,...
				tmp(:,2), [], wshares_kernel_options);

			tmp = 1 - wcumshare_interp.icdf(0.9);
			obj.w_top10share = sfill(tmp, 'w, Top 10% share', 2);

			tmp = 1 - wcumshare_interp.icdf(0.99);
			obj.w_top1share = sfill(tmp, 'w, Top 1% share', 2);

			% Top liquid wealth shares
			values_b = cumsum(obj.grdKFE.b.vec .* obj.pmf_b / obj.liqw.value);
			bcumshare_interp = get_interpolant(values_b,...
				obj.pmf_b, [], wshares_kernel_options);
			tmp = 1 - bcumshare_interp.icdf(0.9);
			obj.lw_top10share = sfill(tmp, 'b, Top 10% share');

			tmp = 1 - bcumshare_interp.icdf(0.99);
			obj.lw_top1share = sfill(tmp, 'b, Top 1% share', 2);

			% Top illiquid wealth shares
			if ~obj.p.OneAsset
				values_a = cumsum(obj.grdKFE.a.vec .* obj.pmf_a(:) / obj.illiqw.value);
				acumshare_interp = get_interpolant(obj.grdKFE.a.vec,...
					obj.pmf_a(:), [], wshares_kernel_options);

				iwshare_interp = @(x) acumshare_interp.icdf(x);
			else
				iwshare_interp = @(x) NaN;
			end

			tmp = 1 - iwshare_interp(0.9);
			obj.iw_top10share = sfill(tmp, 'a, Top 10% share', 2);

			tmp = 1 - iwshare_interp(0.99);
			obj.iw_top1share = sfill(tmp, 'a, Top 1% share', 2);
			
			% Gini coefficient
			tmp = direct_gini(obj.wealth_sorted, obj.pmf_w);
			obj.wgini = sfill(tmp, 'Gini coefficient, wealth');

			% Share of illiquid wealth owned by households in
			% liquid wealth quintiles
			if ~obj.p.OneAsset
				grids = {obj.grdKFE.b.vec, obj.grdKFE.a.vec};
				vals = repmat(shiftdim(obj.grdKFE.a.vec, -1),...
					[obj.p.nb_KFE, 1]);
				integral_a = interp_integral_alt(grids,...
					vals, obj.pmf_b_a);

				amax = obj.grdKFE.a.vec(end);
				lt_b10 = integral_a({obj.lwpercentiles{1}.value, amax});
				lt_b25 = integral_a({obj.lwpercentiles{2}.value, amax});
			else
				lt_b10 = NaN;
				lt_b25 = NaN;
			end

			obj.iwshare_b10 = sfill(lt_b10 / obj.illiqw.value,...
				strcat('Share of illiquid wealth owned by households',...
				' in the bottom 10th percentile of liquid wealth'), 2);
			
			obj.iwshare_b25 = sfill(lt_b25 / obj.illiqw.value,...
				strcat('Share of illiquid wealth owned by households',...
				' in the bottom 25th percentile of liquid wealth'), 2);
		end

		function compute_constrained(obj)
            import HACTLib.aux.multi_sum

            wy_kernel_options = obj.kernel_options;
			wy_kernel_options.x_transform = @(x) log(0.001 + 100*x);
		    wy_kernel_options.x_revert = @(z) (exp(z) - 0.001) / 100;

		    neps = numel(obj.p.epsilon_HtM);
		    obj.constrained_illiq = cell(1, neps);
		    obj.constrained_illiq_pct = cell(1, neps);
		    obj.constrained_illiq_dollars = cell(1, neps);
		    obj.constrained_liq = cell(1, neps);
		    obj.constrained_liq_pct = cell(1, neps);
		    obj.constrained_liq_dollars = cell(1, neps);
		    obj.constrained = cell(1, neps);
		    obj.constrained_pct = cell(1, neps);
		    obj.constrained_dollars = cell(1, neps);
		    for ip = 1:neps
				htm = obj.p.epsilon_HtM(ip);

				tmp = obj.lw_interp.cdf(htm);
				obj.constrained_liq{ip} = sfill(tmp,...
					sprintf('b <= %g', htm));

				obj.constrained_liq_pct{ip} = sfill(tmp,...
					sprintf('b <= %g%% mean ann inc', htm * 100));

				tmp = obj.iw_interp.cdf(htm);
				obj.constrained_illiq{ip} = sfill(tmp,...
					sprintf('a <= %g', htm), 2);

				obj.constrained_illiq_pct{ip} = sfill(tmp,...
					sprintf('a <= %g%% mean ann inc', htm * 100), 2);

				tmp = obj.w_interp.cdf(htm);
				obj.constrained{ip} = sfill(tmp,...
					sprintf('w <= %g', htm), 2);

				obj.constrained_pct{ip} = sfill(tmp,...
					sprintf('w <= %g%% mean ann inc', htm * 100), 2);
			end

			ndollars = numel(obj.p.dollars_HtM);
			for ip = 1:ndollars
				if ~isempty(obj.p.numeraire_in_dollars)
					htm = obj.p.dollars_HtM(ip) / obj.p.numeraire_in_dollars;
					illiq_htm = obj.iw_interp.cdf(htm);
					liq_htm = obj.lw_interp.cdf(htm);
					w_htm = obj.w_interp.cdf(htm);
				else
					illiq_htm = NaN;
					liq_htm = NaN;
					w_htm = NaN;
				end

				obj.constrained_illiq_dollars{ip} = sfill(illiq_htm,...
					sprintf('a <= $%g', obj.p.dollars_HtM(ip)), 2);

				obj.constrained_liq_dollars{ip} = sfill(liq_htm,...
					sprintf('b <= $%g', obj.p.dollars_HtM(ip)));

				obj.constrained_dollars{ip} = sfill(w_htm,...
					sprintf('w <= $%g', obj.p.dollars_HtM(ip)), 2);
			end

			% Fraction paying illiquid asset tax
			z = obj.p.illiquid_tax_threshold;
			if isfinite(z)
				tmp = 1 - obj.iw_interp.cdf(z);
			else
				tmp = 0;
			end
			obj.hhs_paying_wealth_tax = sfill(tmp,...
				'HHs paying tax on illiquid returns', 2);

			% Wealth / (quarterly earnings) < epsilon
			wy_ratio = obj.wealthmat ./ obj.income.y.wide;
			pmf_wy = sum(obj.pmf, 3);

			tmp = sortrows([wy_ratio(:), pmf_wy(:)]);
			values = tmp(:,1);
			keep = values <= 0.5;
			values = values(keep);
			pmf_wy = tmp(keep,2);

			wy_interp = get_interpolant(values, pmf_wy,...
				[], wy_kernel_options);
			
			tmp = wy_interp.cdf(1/6);
			obj.w_lt_ysixth = sfill(...
				tmp, 'w_i <= y_i / 6, biweekly earnings', 2);
			
			tmp = wy_interp.cdf(1/12);
			obj.w_lt_ytwelfth = sfill(...
				tmp, 'w_i <= y_i / 12, weekly earnings', 2);

			% Liquid wealth / (quarterly earnings) < epsilon
			by_ratio = obj.grdKFE.b.vec ./ obj.income.y.wide;
			pmf_by = multi_sum(obj.pmf, [2, 3]);

			tmp = sortrows([by_ratio(:), pmf_by(:)]);
			values = tmp(:,1);
			keep = values <= 0.5;
			values = values(keep);
			pmf_by = tmp(keep,2);

			by_interp = get_interpolant(values, pmf_by,...
				[], wy_kernel_options);
			
			tmp = by_interp.cdf(1/6);
			obj.liqw_lt_ysixth = sfill(...
				tmp, 'b_i <= y_i / 6, biweekly earnings');
            
            tmp = by_interp.cdf(1/12);
			obj.liqw_lt_ytwelfth = sfill(...
				tmp, 'b_i <= y_i / 12, weekly earnings');
            
			% HtM Ratios
			tmp = 1 - obj.w_lt_ysixth.value / obj.liqw_lt_ysixth.value;
			obj.WHtM_over_HtM_biweekly = sfill(tmp,...
				'P(WHtM) / P(HtM), biweekly earnings (y/6)', 2);

			tmp = 1 - obj.w_lt_ytwelfth.value / obj.liqw_lt_ytwelfth.value;
			obj.WHtM_over_HtM_weekly = sfill(tmp,...
				'P(WHtM) / P(HtM), weekly earnings (y/12)', 2);
		end

		function compute_deposit_stats(obj)
			sfill2 = @(y) sfill(NaN, y, 2);

			obj.adjcosts = struct();

			% Adj cost parameters
			obj.adjcosts.kappa0 = sfill2(...
				'kappa0, adj cost coeff on first (linear) term');
			obj.adjcosts.kappa1 = sfill2(...
				'kappa1, adj cost coeff on second (power) term');
			obj.adjcosts.kappa2 = sfill2(...
				'kappa2, power on second term');
			obj.adjcosts.kappa_var = sfill2(...
				'kappa1 ^(-1/kappa2), low values cause mass at high a');
			obj.adjcosts.a_lb = sfill2(...
				'a_lb, parameter s.t. max(a, a_lb) used for adj cost');
			obj.adjcosts.adj_cost_fn = sfill2(...
				'Adjustment cost function');
			obj.adjcosts.mean_cost = sfill2(...
		        'Mean adjustment cost, E[chi(d, a)]');
			obj.adjcosts.mean_d_div_a = sfill2(...
		        'Mean ratio of deposits to assets, E[|d| / max(a, a_lb)]');

			if ~obj.p.OneAsset
				obj.adjcosts.kappa0.value = obj.p.kappa0;
				obj.adjcosts.kappa1.value = obj.p.kappa1;
				obj.adjcosts.kappa2.value = obj.p.kappa2;
				obj.adjcosts.kappa_var.value = obj.p.kappa1 ^ (-1/obj.p.kappa2);
				obj.adjcosts.a_lb.value = obj.p.a_lb;

				lhs = "cost(d,a)";
				term1 = sprintf("%g |d|", obj.p.kappa0);
				term2 = sprintf("(%g / (1 + %g)) |d / max(a,%g)| ^ (1 + %g) * max(a,%g)",...
					obj.p.kappa1, obj.p.kappa2, obj.p.a_lb, obj.p.kappa2, obj.p.a_lb);
				fn_form = strcat(lhs, " = ", term1, " + ", term2);
				obj.adjcosts.adj_cost_fn.value = fn_form;

				% Adj cost statistics
				adj_cost_obj = HACTLib.aux.AdjustmentCost();
				adj_cost_obj.set_from_params(obj.p);
				chii = adj_cost_obj.compute_cost(obj.model.d,...
		        	obj.grdKFE.a.wide);
				obj.adjcosts.mean_cost.value = obj.expectation(chii);

				% Mean abs(d)/a
		        d_div_a = abs(obj.model.d ./ max(obj.grdKFE.a.wide,...
		        	obj.p.a_lb));
		        obj.adjcosts.mean_d_div_a.value = obj.expectation(d_div_a);
			end
		end

		function out = expectation(obj, vals)
			import HACTLib.aux.repmat_auto
			if numel(vals(:)) == numel(obj.pmf(:))
				out = dot(obj.pmf(:), vals(:));
			else
				tmp = repmat_auto(vals, size(obj.pmf));
				out = dot(obj.pmf(:), tmp(:));
			end
		end

		function [pmf_x, cdf_x] = marginal_dists(obj, dims)
			import HACTLib.aux.multi_sum

			sum_dims = 1:4;
			sum_dims = sum_dims(~ismember(sum_dims, dims));

			flatten = true;
			pmf_x = multi_sum(obj.pmf, sum_dims, flatten);

			if isequal(dims, [1, 2])
				dist_sorted = sortrows([obj.wealthmat(:), pmf_x]);
				pmf_x = dist_sorted(:,2);
			end
			cdf_x = cumsum(pmf_x);
		end
	end

	methods (Static)
		function stats = test(mat_path)
			import HACTLib.model_objects.Statistics
			if nargin == 0
				mat_dir = '/home/brian/Documents/temp';
				mat_file = 'stats_vars_2asset.mat'; % stats_vars.mat'; % stats_vars_2asset.mat
				mat_path = fullfile(mat_dir, mat_file);
			end

			load(mat_path);

			stats = Statistics(p, income, grdKFE, KFE);
			stats.compute_statistics();
		end
	end
end

function out = sfill(value, label, asset_indicator)
	% 0 - both assets
	% 1 asset only
	% 2 asset only

	if nargin < 3
		asset_indicator = 0;
	end

	out = struct(...
		'value', value,...
		'label', label,...
		'indicator', asset_indicator...
	);
end

function interp_obj = get_interpolant(values, pmf, dims_to_keep,...
	kernel_options)
	import HACTLib.computation.InterpObj

	already_sorted = isequal(dims_to_keep, []);

	interp_obj = InterpObj();
	interp_obj.set_dist(values, pmf, dims_to_keep, already_sorted);

	interp_obj.configure(kernel_options);
end