function compute_inequality(obj)
	import HACTLib.aux.interp_integral_alt
	import HACTLib.aux.direct_gini
	import HACTLib.aux.multi_sum

	% Top wealth shares
	pmf_w = multi_sum(obj.pmf, [3, 4]);
	tmp = sortrows([obj.wealthmat(:), pmf_w(:)]);
	values_w = cumsum(tmp(:,1) .* tmp(:,2) / obj.totw.value);
	wcumshare_interp = obj.get_interpolant(obj.kernel_options,...
		values_w, tmp(:,2));

	tmp = 1 - wcumshare_interp.icdf(0.9);
	obj.w_top10share = obj.sfill(tmp, 'w, Top 10% share', 2);

	tmp = 1 - wcumshare_interp.icdf(0.99);
	obj.w_top1share = obj.sfill(tmp, 'w, Top 1% share', 2);

	% Top liquid wealth shares
	values_b = cumsum(obj.grdKFE.b.vec .* obj.pmf_b / obj.liqw.value);
	bcumshare_interp = obj.get_interpolant(obj.kernel_options,...
		values_b, obj.pmf_b);
	tmp = 1 - bcumshare_interp.icdf(0.9);
	obj.lw_top10share = obj.sfill(tmp, 'b, Top 10% share');

	tmp = 1 - bcumshare_interp.icdf(0.99);
	obj.lw_top1share = obj.sfill(tmp, 'b, Top 1% share', 2);

	% Top illiquid wealth shares
	if ~obj.p.OneAsset
		values_a = cumsum(obj.grdKFE.a.vec .* obj.pmf_a(:) / obj.illiqw.value);
		acumshare_interp = obj.get_interpolant(obj.kernel_options,...
			obj.grdKFE.a.vec, obj.pmf_a(:));

		iwshare_interp = @(x) acumshare_interp.icdf(x);
	else
		iwshare_interp = @(x) NaN;
	end

	tmp = 1 - iwshare_interp(0.9);
	obj.iw_top10share = obj.sfill(tmp, 'a, Top 10% share', 2);

	tmp = 1 - iwshare_interp(0.99);
	obj.iw_top1share = obj.sfill(tmp, 'a, Top 1% share', 2);
	
	% Gini coefficient
	tmp = direct_gini(obj.wealth_sorted, obj.pmf_w);
	obj.wgini = obj.sfill(tmp, 'Gini coefficient, wealth');

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

	obj.iwshare_b10 = obj.sfill(lt_b10 / obj.illiqw.value,...
		strcat('Share of illiquid wealth owned by households',...
		' in the bottom 10th percentile of liquid wealth'), 2);
	
	obj.iwshare_b25 = obj.sfill(lt_b25 / obj.illiqw.value,...
		strcat('Share of illiquid wealth owned by households',...
		' in the bottom 25th percentile of liquid wealth'), 2);
end