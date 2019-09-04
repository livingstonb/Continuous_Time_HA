function [H_special,c_special,d_special] = deal_with_special_case(p,income,grids,r_b_mat,VaB)
	% this function computes the policy functions associated
	% with the special case that b = bmin, a > 0 and households
	% withdraw only enough to consume, so that a does not accumulate
	
	lambda = zeros(p.nb,p.na,p.ny,p.nz);
	rb = r_b_mat(1,1,1,1) * p.bmin;
	for ia = 2:p.na
	for k = 1:p.ny
	for iz = 1:p.nz
		l_fn = @(x) aux.lambda_function(x,VaB(1,ia,k,iz),rb,grids.a.vec(ia),income.y.vec(k),p);

		options = optimset('TolX',1e-8,'TolFun',1e-11);
        [lambda(1,ia,k,iz),fval] = fminbnd(@(x) l_fn(x)^2,0,1e9,options);

        if fval > 1e-10
            error('no convergence for special case')
        end
	end
	end
	end

	chi1inv_arg = VaB ./ lambda - 1;
    d_special = aux.adj_cost_deriv_inverse(chi1inv_arg,grids.a.matrix,p);
    d_special(~isfinite(d_special)) = 0;


    c_special = zeros(p.nb,p.na,p.ny,p.nz);
    c_special(1,2:p.na,:,:) = lambda(1,2:p.na,:,:) .^(-1/p.riskaver);

    H_special = - 1e12 * ones(p.nb,p.na,p.ny,p.nz);
    H_special(1,2:p.na,:,:) = VaB(1,2:p.na,:,:) .* d_special(1,2:p.na,:,:);
end