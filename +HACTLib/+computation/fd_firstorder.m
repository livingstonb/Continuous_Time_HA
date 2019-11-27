function [fdB, fdF] = fd_firstorder(values, deltaB, deltaF, dim)
	% Approximates the derivative of the first argument
	% via a first-order finite difference. Zeros are
	% placed at the end of the dimension 'dim' for the
	% forward difference and at the beginning of 'dim'
	% for the backward difference.

	zeros_shape = size(values);
	zeros_shape(dim) = 1;

	diffs = diff(values, 1, dim);
	fdB = cat(dim, zeros(zeros_shape), diffs) ./ deltaB;
	fdF = cat(dim, diffs, zeros(zeros_shape)) ./ deltaF;
end