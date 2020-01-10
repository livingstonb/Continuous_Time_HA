function check_for_required_properties(object, props)
	% Checks 'object' to verify that it has at least
	% the required attributes.
	
	assert(iscell(props), "Second arg to required attributes must be cell array")

	for i = 1:numel(props)
		if ~isprop(object, props{i})
			msg = sprintf("Required attribute '%s' not found", props{i});
			error(msg)
		end
	end
end