function fd = fd_firstorder(values, deltaF, deltaB, dim)
	% Approximates the derivative of the first argument
	% via a first-order finite difference. Zeros are
	% placed at the end of the dimension 'dim' for the
	% forward difference and at the beginning of 'dim'
	% for the backward difference.

	zeros_shape = size(values);
	zeros_shape(dim) = 1;

	diffs = diff(values, 1, dim);	
	fd.F = cat(dim, diffs, zeros(zeros_shape)) ./ deltaF;
	fd.B = cat(dim, zeros(zeros_shape), diffs) ./ deltaB;
end