function assert_correct_shape(object, shape)
	% Performs assertions that an input is an
	% array with shape equal to 'shape'.

	assert(isnumeric(object) || iscell(object),...
		"Input is not an array")

	assert(isequal(size(object), shape),...
		"Input is not of the declared shape")
end