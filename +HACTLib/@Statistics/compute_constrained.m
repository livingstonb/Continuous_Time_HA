function compute_constrained(obj)
	import HACTLib.aux.multi_sum

	lw_interp = griddedInterpolant(obj.bgrid, cumsum(obj.pmf_b), 'pchip', 'nearest');
	iw_interp = griddedInterpolant(obj.agrid, cumsum(obj.pmf_a), 'pchip', 'nearest');

	sorted_mat = sortrows([obj.wealthmat(:), obj.pmf_w(:)]);
	[w_u, iu] = unique(sorted_mat(:,1), 'last');
	cdf_w = cumsum(sorted_mat(:,2));
	w_interp = griddedInterpolant(w_u, cdf_w(iu), 'pchip', 'nearest');

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

		tmp = lw_interp(htm);
		obj.constrained_liq{ip} = obj.sfill(tmp,...
			sprintf('b <= %g', htm));

		obj.constrained_liq_pct{ip} = obj.sfill(tmp,...
			sprintf('b <= %g%% mean ann inc', htm * 100));

		tmp = iw_interp(htm);
		obj.constrained_illiq{ip} = obj.sfill(tmp,...
			sprintf('a <= %g', htm), 2);

		obj.constrained_illiq_pct{ip} = obj.sfill(tmp,...
			sprintf('a <= %g%% mean ann inc', htm * 100), 2);

		tmp = w_interp(htm);
		obj.constrained{ip} = obj.sfill(tmp,...
			sprintf('w <= %g', htm), 2);

		obj.constrained_pct{ip} = obj.sfill(tmp,...
			sprintf('w <= %g%% mean ann inc', htm * 100), 2);
	end

	ndollars = numel(obj.p.dollars_HtM);
	for ip = 1:ndollars
		if ~isempty(obj.p.numeraire_in_dollars)
			htm = obj.p.dollars_HtM(ip) / obj.p.numeraire_in_dollars;
			illiq_htm = iw_interp(htm);
			liq_htm = lw_interp(htm);
			w_htm = w_interp(htm);
		else
			illiq_htm = NaN;
			liq_htm = NaN;
			w_htm = NaN;
		end

		obj.constrained_illiq_dollars{ip} = obj.sfill(illiq_htm,...
			sprintf('a <= $%g', obj.p.dollars_HtM(ip)), 2);

		obj.constrained_liq_dollars{ip} = obj.sfill(liq_htm,...
			sprintf('b <= $%g', obj.p.dollars_HtM(ip)));

		obj.constrained_dollars{ip} = obj.sfill(w_htm,...
			sprintf('w <= $%g', obj.p.dollars_HtM(ip)), 2);
	end

	% Liquid wealth / (quarterly earnings) < epsilon
	kernel_options = struct();
	kernel_options.ktype = 'triweight';
    kernel_options.h = 0.3;
    kernel_options.rescale_and_log = true;
    kernel_options.force_fit_cdf_low = [];

	by_ratio = obj.bgrid ./ obj.income.y.wide;
	pmf_by = multi_sum(obj.pmf, [2, 3]);

	tmp = sortrows([by_ratio(:), pmf_by(:)]);
	by_interp = get_interpolant(kernel_options,...
		tmp(:,1), tmp(:,2), 0.3, []);
	
	obj.liqw_lt_ysixth = obj.sfill(by_interp.cdf(1/6), 'b_i <= y_i / 6');
   	obj.liqw_lt_ytwelfth = obj.sfill(by_interp.cdf(1/12), 'b_i <= y_i / 12');

	% Wealth / (quarterly earnings) < epsilon
	wy_ratio = obj.wealthmat ./ obj.income.y.wide;
	pmf_wy = sum(obj.pmf, 3);

	tmp = sortrows([wy_ratio(:), pmf_wy(:)]);
	wy_interp = get_interpolant(kernel_options,...
		tmp(:,1), tmp(:,2), 0.3);
	
	obj.w_lt_ysixth = obj.sfill(wy_interp.cdf(1/6), 'w_i <= y_i / 6', 2);
	obj.w_lt_ytwelfth = obj.sfill(wy_interp.cdf(1/12), 'w_i <= y_i / 12', 2);

	% HtM Ratios
	tmp = 1 - obj.w_lt_ysixth.value / obj.liqw_lt_ysixth.value;
	obj.WHtM_over_HtM_biweekly = obj.sfill(tmp,...
		'P(WHtM) / P(HtM), HtM in terms of y/6', 2);

	tmp = 1 - obj.w_lt_ytwelfth.value / obj.liqw_lt_ytwelfth.value;
	obj.WHtM_over_HtM_weekly = obj.sfill(tmp,...
		'P(WHtM) / P(HtM), HtM in terms of y/12', 2);
end

function interp_obj = get_interpolant(kernel_options, varargin)
	import HACTLib.computation.InterpObj
	interp_obj = InterpObj();
	interp_obj.set_dist(varargin{:});
	interp_obj.configure(kernel_options);
end