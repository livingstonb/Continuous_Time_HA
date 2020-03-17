function add_mpcs(obj, mpc_obj, simulated)
	if mpc_obj.options.liquid_mpc
		asset_indicator = 0;
		mpc_type = '';
	else
		asset_indicator = 2;
		mpc_type = 'illiq';
	end
	sfill2 = @(x,y) obj.sfill(x, y, asset_indicator);

	empty_stat = sfill2([], []);

	empty_mpc_struct = struct(...
		'shock_normalized', empty_stat,...
		'shock', empty_stat,...
		'quarterly', empty_stat,...
		'annual', empty_stat...
	);

	nshocks = numel(obj.p.mpc_shocks);
	for ishock = 1:nshocks
		shock = obj.p.mpc_shocks(ishock);
		shock_label = obj.p.quantity2label(shock);
		mpcs_stats(ishock) = empty_mpc_struct;

		mpcs_stats(ishock).shock_normalized = sfill2(...
			shock * 100, 'Shock size, (% of mean ann inc)');

		mpcs_stats(ishock).shock = sfill2(shock_label,...
			'Shock size');

		tmp = 100 * mpc_obj.mpcs(ishock).quarterly(1);
		label = sprintf(...
			'Quarterly %s MPC (%%), out of %s',...
			mpc_type, shock_label);
		mpcs_stats(ishock).quarterly = sfill2(tmp, label);

		tmp = 100 * mpc_obj.mpcs(ishock).annual;
		label = sprintf(...
			'Annual %s MPC (%%), out of %s',...
			mpc_type, shock_label);
		mpcs_stats(ishock).annual = sfill2(tmp, label);
	end

	if mpc_obj.options.liquid_mpc
		obj.mpcs = mpcs_stats;

		obj.mpcs_over_ss = cell(1, nshocks);
		for ishock = 1:nshocks
			obj.mpcs_over_ss{ishock} =  mpc_obj.mpcs(ishock).mpcs(:,1);
		end
	else
		obj.illiquid_mpcs = mpcs_stats;
	end
end