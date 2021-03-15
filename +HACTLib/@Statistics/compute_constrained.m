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
    kopts = obj.kernel_options;
    ktype0 = kopts.ktype;
    h0 = kopts.h;
    rescale_and_log0 = kopts.rescale_and_log;

    kopts.ktype = 'gaussian';
    kopts.h = 0.15;
%     kopts.force_fit_cdf_low = [0.07, 0.5];
    kopts.rescale_and_log = true;

	by_ratio = obj.bgrid ./ obj.income.y.wide;
	pmf_by = multi_sum(obj.pmf, [2, 3]);

	tmp = sortrows([by_ratio(:), pmf_by(:)]);
	by_interp = obj.get_interpolant(kopts,...
		tmp(:,1), tmp(:,2), 0.4, []);
	
	tmp = by_interp.cdf(1/6);
	obj.liqw_lt_ysixth = obj.sfill(...
		tmp, 'b_i <= y_i / 6');
    
    tmp = by_interp.cdf(1/12);
	obj.liqw_lt_ytwelfth = obj.sfill(...
		tmp, 'b_i <= y_i / 12');

	% Wealth / (quarterly earnings) < epsilon
	wy_ratio = obj.wealthmat ./ obj.income.y.wide;
	pmf_wy = sum(obj.pmf, 3);

	tmp = sortrows([wy_ratio(:), pmf_wy(:)]);
	wy_interp = obj.get_interpolant(kopts,...
		tmp(:,1), tmp(:,2), 0.4);
	
	tmp = wy_interp.cdf(1/6);
	obj.w_lt_ysixth = obj.sfill(...
		tmp, 'w_i <= y_i / 6', 2);
	
	tmp = wy_interp.cdf(1/12);
	obj.w_lt_ytwelfth = obj.sfill(...
		tmp, 'w_i <= y_i / 12', 2);

	% HtM Ratios
	tmp = 1 - obj.w_lt_ysixth.value / obj.liqw_lt_ysixth.value;
	obj.WHtM_over_HtM_biweekly = obj.sfill(tmp,...
		'P(WHtM) / P(HtM), HtM in terms of y/6', 2);

	tmp = 1 - obj.w_lt_ytwelfth.value / obj.liqw_lt_ytwelfth.value;
	obj.WHtM_over_HtM_weekly = obj.sfill(tmp,...
		'P(WHtM) / P(HtM), HtM in terms of y/12', 2);
    
    kopts.rescale_and_log = rescale_and_log0;
    kopts.ktype = ktype0;
    kopts.h = h0;
end