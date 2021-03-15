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
	obj.adjcosts.mean_chi_div_d = obj.sfill(NaN,...
		'E[chi(d,a)/abs(d) | d != 0]', 2);

	npct = numel(obj.p.wpercentiles);
	obj.adjcosts.chi_div_d_pctiles = cell(1, npct);
	for ip = 1:npct
		pct_at = obj.p.wpercentiles(ip);
		obj.adjcosts.chi_div_d_pctiles{ip} = obj.sfill(NaN,...
			sprintf('chi/abs(d), %gth pctile condl on d != 0', pct_at), 2);
	end

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
        	shiftdim(obj.agrid, -1));
		obj.adjcosts.mean_cost.value = obj.expectation(chii);

		% Mean abs(d) / a
        d_div_a = abs(obj.model.d ./ max(...
        	shiftdim(obj.agrid, -1), obj.p.a_lb));
        obj.adjcosts.mean_d_div_a.value = obj.expectation(d_div_a);

        % Conditional distribution of chi / |d|
        chii = adj_cost_obj.compute_cost(obj.model.d,...
        	shiftdim(obj.agrid, -1));
        chi_div_d = chii ./ abs(obj.model.d);
        nonzdeposits = abs(obj.model.d) > 1.0e-7;
        chi_div_d = chi_div_d(:);
        chi_div_d = chi_div_d(nonzdeposits(:));
        condpmf = obj.pmf(:);
        condpmf = condpmf(nonzdeposits(:));
        condpmf = condpmf / sum(condpmf);

        interpolant = HACTLib.aux.pctile_interpolant(chi_div_d, condpmf);
        obj.adjcosts.mean_chi_div_d.value = condpmf' * chi_div_d;
		for ip = 1:npct
			pct_at = obj.p.wpercentiles(ip);
			obj.adjcosts.chi_div_d_pctiles{ip}.value = interpolant(pct_at / 100);
		end
	end
end