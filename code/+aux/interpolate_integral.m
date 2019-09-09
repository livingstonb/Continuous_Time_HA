function interpolant = interpolate_integral(gridValues, integrandValues, pmf)
	% returns an interpolant that interpolates to find the value of
	% int_0^{epsilon} values(a)g(a)da for a given epsilon
	%
	% evaluate the above integral by calling interpolant(epsilon)

	sortedInputs = sortrows([gridValues(:) integrandValues(:) pmf(:)]);
	gridSorted = sortedInputs(:,1);
	integrandSorted = sortedInputs(:,2);
	pmfSorted = sortedInputs(:,3);

	integralValues = cumsum(integrandSorted .* pmfSorted);

	[gridUnique,uniqueInds] = unique(gridSorted,'last');
	integralUnique = integralValues(uniqueInds);

	interpolant = griddedInterpolant(gridUnique,integralUnique,'linear');

end