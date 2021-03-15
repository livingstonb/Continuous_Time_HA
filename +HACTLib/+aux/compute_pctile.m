function pctile = compute_pctile(values, pmf)
	sorted_mat = sortrows([values(:), pmf_x(:)]);
	[cdf_u, iu] = unique(cumsum(sorted_mat(:,2)), 'last');
	values_u = sorted_mat(iu,1);

	pctile = interp1d(cdf_u, values_u, 'pchip', 'nearest');
end