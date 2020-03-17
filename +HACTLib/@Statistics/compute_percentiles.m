function compute_percentiles(obj)
	npct = numel(obj.p.wpercentiles);
    obj.lwpercentiles = cell(1, npct);
    obj.iwpercentiles = cell(1, npct);
    obj.wpercentiles = cell(1, npct);
	for ip = 1:npct
		pct_at = obj.p.wpercentiles(ip);

		tmp_b = sprintf('b, %gth pctile', pct_at);
		obj.lwpercentiles{ip} = obj.sfill(...
			obj.lw_interp.icdf(pct_at/100), tmp_b);

		tmp_a = sprintf('a, %gth pctile', pct_at);
		obj.iwpercentiles{ip} = obj.sfill(...
			obj.iw_interp.icdf(pct_at/100), tmp_a, 2);

		tmp_w = sprintf('w, %gth pctile', pct_at);
		obj.wpercentiles{ip} = obj.sfill(...
			obj.w_interp.icdf(pct_at/100), tmp_w, 2);
	end

	obj.median_liqw = obj.sfill(obj.lw_interp.icdf(0.5),...
		'b, median');
	obj.median_illiqw = obj.sfill(obj.iw_interp.icdf(0.5),...
		'a, median', 2);
	obj.median_totw = obj.sfill(obj.w_interp.icdf(0.5),...
		'w, median', 2);
end