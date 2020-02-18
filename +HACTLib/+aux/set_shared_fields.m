function new_params_struct = set_shared_fields(params, shared_params)
	pfields = fields(shared_params)';

	
	for ii = 1:numel(params)
		new_params{ii} = params(ii);
		for pfield = pfields
			val = shared_params.(pfield{1});
			new_params{ii} = setfield(new_params{ii}, pfield{1}, val);
		end
	end

	new_params_struct = new_params{1};
	for ii = 2:numel(params)
		new_params_struct(ii) = new_params{ii};
	end
end