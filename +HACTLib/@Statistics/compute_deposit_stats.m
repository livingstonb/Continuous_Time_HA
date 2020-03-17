function compute_deposit_stats(obj)
	sfill2 = @(y) obj.sfill(NaN, y, 2);

	obj.adjcosts = struct();

	% Adj cost parameters
	obj.adjcosts.kappa0 = sfill2(...
		'kappa0, adj cost coeff on first (linear) term');
	obj.adjcosts.kappa1 = sfill2(...
		'kappa1, adj cost coeff on second (power) term');
	obj.adjcosts.kappa2 = sfill2(...
		'kappa2, power on second term');
	obj.adjcosts.kappa_var = sfill2(...
		'kappa1 ^(-1/kappa2), low values cause mass at high a');
	obj.adjcosts.a_lb = sfill2(...
		'a_lb, parameter s.t. max(a, a_lb) used for adj cost');
	obj.adjcosts.adj_cost_fn = sfill2(...
		'Adjustment cost function');
	obj.adjcosts.mean_cost = sfill2(...
        'Mean adjustment cost, E[chi(d, a)]');
	obj.adjcosts.mean_d_div_a = sfill2(...
        'Mean ratio of deposits to assets, E[|d| / max(a, a_lb)]');

	if ~obj.p.OneAsset
		obj.adjcosts.kappa0.value = obj.p.kappa0;
		obj.adjcosts.kappa1.value = obj.p.kappa1;
		obj.adjcosts.kappa2.value = obj.p.kappa2;
		obj.adjcosts.kappa_var.value = obj.p.kappa1 ^ (-1/obj.p.kappa2);
		obj.adjcosts.a_lb.value = obj.p.a_lb;

		lhs = "cost(d,a)";
		term1 = sprintf("%g |d|", obj.p.kappa0);
		term2 = sprintf("(%g / (1 + %g)) |d / max(a,%g)| ^ (1 + %g) * max(a,%g)",...
			obj.p.kappa1, obj.p.kappa2, obj.p.a_lb, obj.p.kappa2, obj.p.a_lb);
		fn_form = strcat(lhs, " = ", term1, " + ", term2);
		obj.adjcosts.adj_cost_fn.value = fn_form;

		% Adj cost statistics
		adj_cost_obj = HACTLib.aux.AdjustmentCost();
		adj_cost_obj.set_from_params(obj.p);
		chii = adj_cost_obj.compute_cost(obj.model.d,...
        	obj.grdKFE.a.wide);
		obj.adjcosts.mean_cost.value = obj.expectation(chii);

		% Mean abs(d)/a
        d_div_a = abs(obj.model.d ./ max(obj.grdKFE.a.wide,...
        	obj.p.a_lb));
        obj.adjcosts.mean_d_div_a.value = obj.expectation(d_div_a);
	end
end