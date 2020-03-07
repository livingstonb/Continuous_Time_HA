function fraction_below = find_constrained(values, pmf, thresholds)
	% given values and a pmf, this functions finds the fraction
	% of households below given thresholds
	sorted_by_values = sortrows([values(:) pmf(:)]);

	values_sorted = sorted_by_values(:,1);
	cdf_sorted = cumsum(sorted_by_values(:,2));

	[values_unique, iu] = unique(values_sorted, 'last');
	cdf_unique = cdf_sorted(iu);

	values_interp = griddedInterpolant(values_unique, cdf_unique,...
		'pchip', 'nearest');
	fraction_below = values_interp(thresholds);
end