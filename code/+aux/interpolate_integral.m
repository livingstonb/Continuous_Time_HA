function interpolant = interpolate_integral(gridValues, integrandValues, pmf)
	% creates an interpolant that approximates the value of the integral
	% int_0^{epsilon} values(a)g(a)da for a given epsilon

	% Parameters
	% ----------
	% gridValues : values at which the integrand is evaluated
	%
	% integrandValues : values of the integrand
	%
	% pmf : the probability mass function over states
	%
	% Results
	% -------
	% interpolant : a griddedInterpolant object such that interpolant(x)
	%	is the approximated value of the integral from 0 to x

	sortedInputs = sortrows([gridValues(:) integrandValues(:) pmf(:)]);
	gridSorted = sortedInputs(:,1);
	integrandSorted = sortedInputs(:,2);
	pmfSorted = sortedInputs(:,3);

	integralValues = cumsum(integrandSorted .* pmfSorted);

	[gridUnique,uniqueInds] = unique(gridSorted,'last');
	integralUnique = integralValues(uniqueInds);

	interpolant = griddedInterpolant(gridUnique,integralUnique,'linear');

end