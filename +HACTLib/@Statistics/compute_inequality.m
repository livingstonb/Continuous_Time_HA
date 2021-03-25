function compute_inequality(obj)
	import HACTLib.aux.interp_integral_alt
	import HACTLib.aux.direct_gini
	import HACTLib.aux.multi_sum

	% Top wealth shares
	pmf_w = multi_sum(obj.pmf, [3, 4]);
	tmp = sortrows([obj.wealthmat(:), pmf_w(:)]);
	values_w = cumsum(tmp(:,1) .* tmp(:,2) / obj.totw.value);

	[cdf_u, iu] = unique(cumsum(tmp(:,2)), 'last');
	wcumshare_interp = griddedInterpolant(cdf_u, values_w(iu), 'pchip', 'nearest');

	tmp = 1 - wcumshare_interp(0.9);
	obj.w_top10share = obj.sfill(tmp, 'w, Top 10% share', 2);

	tmp = 1 - wcumshare_interp(0.99);
	obj.w_top1share = obj.sfill(tmp, 'w, Top 1% share', 2);

	% Top liquid wealth shares
	values_b = cumsum(obj.bgrid .* obj.pmf_b / obj.liqw.value);
	[cdf_b_u, iu_b] = unique(cumsum(obj.pmf_b), 'last');
	bcumshare_interp = griddedInterpolant(cdf_b_u, values_b(iu_b), 'pchip', 'nearest');

	tmp = 1 - bcumshare_interp(0.9);
	obj.lw_top10share = obj.sfill(tmp, 'b, Top 10% share');

	tmp = 1 - bcumshare_interp(0.99);
	obj.lw_top1share = obj.sfill(tmp, 'b, Top 1% share', 2);

	% Top illiquid wealth shares
	if ~obj.p.OneAsset
		values_a = cumsum(obj.agrid .* obj.pmf_a(:) / obj.illiqw.value);
		[cdf_a_u, iu_a] = unique(cumsum(obj.pmf_a), 'last');
		acumshare_interp = griddedInterpolant(cdf_a_u, values_a(iu_a), 'pchip', 'nearest');

		iwshare_interp = @(x) acumshare_interp(x);
	else
		iwshare_interp = @(x) NaN;
	end

	tmp = 1 - iwshare_interp(0.9);
	obj.iw_top10share = obj.sfill(tmp, 'a, Top 10% share', 2);

	tmp = 1 - iwshare_interp(0.99);
	obj.iw_top1share = obj.sfill(tmp, 'a, Top 1% share', 2);

	% Gini coefficient
	tmp = direct_gini(obj.wealthmat, obj.pmf_w);
	obj.wgini = obj.sfill(tmp, 'Gini coefficient, wealth');

	% Share of illiquid wealth owned by households in
	% liquid wealth quintiles
	if ~obj.p.OneAsset
		grids = {obj.bgrid, obj.agrid};
		vals = repmat(shiftdim(obj.agrid, -1),...
			[obj.p.nb_KFE, 1]);
		integral_a = interp_integral_alt(grids,...
			vals, obj.pmf_b_a);

		amax = obj.agrid(end);
		lt_b10 = integral_a({obj.lwpercentiles{1}.value, amax});
		lt_b25 = integral_a({obj.lwpercentiles{2}.value, amax});
	else
		lt_b10 = NaN;
		lt_b25 = NaN;
	end

	obj.iwshare_b10 = obj.sfill(lt_b10 / obj.illiqw.value,...
		strcat('Share of illiquid wealth owned by households',...
		' in the bottom 10th percentile of liquid wealth'), 2);
	
	obj.iwshare_b25 = obj.sfill(lt_b25 / obj.illiqw.value,...
		strcat('Share of illiquid wealth owned by households',...
		' in the bottom 25th percentile of liquid wealth'), 2);
end