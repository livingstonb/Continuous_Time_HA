function check_for_required_properties(object, props)
	
	assert(iscell(props), "Second arg to required properties must be cell array")

	for i = 1:numel(props)
		if ~isprop(object, props{i})
			msg = sprintf("Required attribute '%s' not found", props{i});
			error(msg)
		end
	end
end