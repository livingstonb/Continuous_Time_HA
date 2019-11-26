function assert_same_shape(varargin)
	% Performs assertions that the inputs are arrays
	% of the same shape.

	if nargin < 2
		return
	end

	for k = 1:nargin
		assert(isnumeric(varargin{k}) || iscell(varargin{k}),...
			"One or more of the inputs is not an array")
	end

	shape = size(varargin{1});
	for k = 2:nargin
		assert(isequal(shape, size(varargin{k})),...
			"Inputs do not have the same shape")
	end
end