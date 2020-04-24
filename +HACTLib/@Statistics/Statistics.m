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

		bgrid;
		agrid;

		nb;
		na;
		nz;
		ny;

		params = struct();

		other = struct();
	end

	properties (Access=protected)
		p;
		income;
		model;

		wealth_sorted;
		wealthmat;

		pmf_w;
        
        kernel_options;

        lw_interp;
        iw_interp;
        w_interp;
	end

	methods
		function obj = Statistics(p, income, grdKFE, model)
			obj.p = p;
			obj.income = income;
			obj.model = model;

			obj.bgrid = grdKFE.b.vec;
			obj.agrid = grdKFE.a.vec;

			obj.na = p.na_KFE;
			obj.nb = p.nb_KFE;
			obj.nz = p.nz;
			obj.ny = income.ny;

			tmp = obj.bgrid + shiftdim(obj.agrid, -1);
			obj.wealthmat = tmp;
			obj.wealth_sorted = sortrows(tmp(:));

			obj.pmf = model.g .* grdKFE.trapezoidal.matrix;
		end

		function compute_statistics(obj, kernel_options)
			if nargin == 2
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
			obj.params.bmax = obj.sfill(obj.p.bmax,...
				'Max liquid assets, parameter');
			obj.params.amax = obj.sfill(obj.p.amax,...
				'Max illiquid assets, parameter', 2);
			obj.params.bequests = obj.sfill(obj.p.Bequests,...
				'Bequests, on or off');
			obj.params.deathrate = obj.sfill(obj.p.deathrate,...
				'Death rate (quarterly)');
			obj.params.r_b = obj.sfill(obj.p.r_b,...
				'Liquid asset return (quarterly)');
			obj.params.r_a = obj.sfill(obj.p.r_a,...
				'Illiquid asset return (quarterly)', 2);
			obj.params.borrowlim = obj.sfill(obj.p.bmin,...
				'Borrowing limit');
			obj.params.riskaver = obj.sfill(obj.p.riskaver,...
				'CRRA coefficient');
			obj.params.numeraire = obj.sfill(obj.p.numeraire_in_dollars,...
				'Value of the numeraire, mean annual earning, in $');
		end

		function clean(obj)
			obj.p = [];
			obj.income = [];
			obj.model = [];
			obj.wealthmat = [];
			obj.wealth_sorted = [];
			obj.pmf_w = [];
			obj.lw_interp = [];
			obj.iw_interp = [];
			obj.w_interp = [];
		end
	end

	methods (Access=protected)
		function compute_intro_stats(obj)
			obj.rho = obj.sfill(obj.p.rho, 'rho');
		    obj.beta_Q = obj.sfill(exp(-obj.p.rho), 'beta (quarterly)');
		    obj.beta_A = obj.sfill(exp(-4 * obj.p.rho), 'beta (annualized)');

		    tmp = obj.expectation(shiftdim(obj.agrid, -1));
		    obj.illiqw = obj.sfill(tmp, 'Mean illiquid wealth', 2);

		    tmp = obj.expectation(obj.bgrid);
		    obj.liqw = obj.sfill(tmp, 'Mean liquid wealth');

		    tmp = obj.expectation(...
		    	obj.bgrid + shiftdim(obj.agrid, -1));
		    obj.totw = obj.sfill(tmp, 'Mean total wealth', 2);

		    tmp = obj.expectation(obj.model.s==0);
		    obj.sav0 = obj.sfill(tmp, 's = 0');
		end

		function construct_distributions(obj)
			import HACTLib.aux.multi_sum

		    [obj.pmf_b, obj.cdf_b] = obj.marginal_dists(1);
		    [obj.pmf_a, obj.cdf_a] = obj.marginal_dists(2);
		    obj.pmf_w = obj.marginal_dists([1, 2]);
		    obj.pmf_b_a = multi_sum(obj.pmf, [3, 4]);

		    obj.lw_interp = obj.get_interpolant(obj.kernel_options,...
		    	obj.bgrid, obj.pmf, [1]);
		    obj.iw_interp = obj.get_interpolant(obj.kernel_options,...
		    	obj.agrid, obj.pmf, [2]);
		    obj.w_interp = obj.get_interpolant(obj.kernel_options,...
		    	obj.wealthmat, obj.pmf, [1, 2]);
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

		function interp_obj = get_interpolant(kernel_options, varargin)
			import HACTLib.computation.InterpObj

			interp_obj = InterpObj();
			interp_obj.set_dist(varargin{:});
			interp_obj.configure(kernel_options);
		end
	end
end