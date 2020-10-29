function compute_constrained(obj)
	import HACTLib.aux.multi_sum

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
		obj.constrained_liq{ip} = obj.sfill(tmp,...
			sprintf('b <= %g', htm));

		obj.constrained_liq_pct{ip} = obj.sfill(tmp,...
			sprintf('b <= %g%% mean ann inc', htm * 100));

		tmp = obj.iw_interp.cdf(htm);
		obj.constrained_illiq{ip} = obj.sfill(tmp,...
			sprintf('a <= %g', htm), 2);

		obj.constrained_illiq_pct{ip} = obj.sfill(tmp,...
			sprintf('a <= %g%% mean ann inc', htm * 100), 2);

		tmp = obj.w_interp.cdf(htm);
		obj.constrained{ip} = obj.sfill(tmp,...
			sprintf('w <= %g', htm), 2);

		obj.constrained_pct{ip} = obj.sfill(tmp,...
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

		obj.constrained_illiq_dollars{ip} = obj.sfill(illiq_htm,...
			sprintf('a <= $%g', obj.p.dollars_HtM(ip)), 2);

		obj.constrained_liq_dollars{ip} = obj.sfill(liq_htm,...
			sprintf('b <= $%g', obj.p.dollars_HtM(ip)));

		obj.constrained_dollars{ip} = obj.sfill(w_htm,...
			sprintf('w <= $%g', obj.p.dollars_HtM(ip)), 2);
	end

	% Fraction paying illiquid asset tax
	z = obj.p.illiquid_tax_threshold;
	if isfinite(z)
		tmp = 1 - obj.iw_interp.cdf(z);
	else
		tmp = 0;
	end
	obj.hhs_paying_wealth_tax = obj.sfill(tmp,...
		'HHs paying tax on illiquid returns', 2);

	% Liquid wealth / (quarterly earnings) < epsilon
    kopts = obj.kernel_options;
    kopts.h = 0.15;

	by_ratio = obj.bgrid ./ obj.income.y.wide;
	pmf_by = multi_sum(obj.pmf, [2, 3]);

	tmp = sortrows([by_ratio(:), pmf_by(:)]);
	by_interp = obj.get_interpolant(kopts,...
		tmp(:,1), tmp(:,2));
	
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
		tmp(:,1), tmp(:,2));
	
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
end