function interpolant = interpolate_integral(grid_values, integrand_values, pmf)
	% Creates an interpolant that approximates the value of the integral
	% int_0^{epsilon} values(a)g(a)da for a given epsilon.
	%
	% Parameters
	% ----------
	% gridValues : Values at which the integrand is evaluated.
	%
	% integrand_values : Values of the integrand.
	%
	% pmf : The probability mass function over states.
	%
	% Results
	% -------
	% interpolant : A griddedInterpolant object such that interpolant(x)
	%	is the approximated value of the integral from 0 to x.

	validate_inputs(grid_values, integrand_values, pmf);

	sorted_inputs = sortrows([grid_values(:) integrand_values(:) pmf(:)]);
	grid_sorted = sorted_inputs(:,1);
	integrand_sorted = sorted_inputs(:,2);
	pmf_sorted = sorted_inputs(:,3);

	integral_values = cumsum(integrand_sorted .* pmf_sorted);

	[grid_unique, unique_inds] = unique(grid_sorted,'last');
	integral_unique = integral_values(unique_inds);

	interpolant = griddedInterpolant(...
		grid_unique, integral_unique, 'linear');
end

function validate_inputs(grid_values, integrand_values, pmf)
	
end