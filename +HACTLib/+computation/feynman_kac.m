function cumcon_update = feynman_kac(p, grids, income, cumcon_t, FKmats, c, stepsize)
	% Computes inflows for death
	%
	% Parameters
	% ----------
	% cumcon_t_k : cumulative consumption for period t,
	%	of shape (nb_KFE*na_KFE*nz, ny)
	
	cumcon_update = zeros(size(cumcon_t));
	cumcon_t_k = reshape(cumcon_t, [], income.ny);
	for k = 1:income.ny
		indx_k = ~ismember(1:income.ny, k);
		ytrans_cc_k = sum(income.ytrans(k,indx_k) .* cumcon_t_k(:,indx_k), 2);

        if p.Bequests
        	deathin_cc_k = p.deathrate * cumcon_t_k(:,k);
        else
        	cumcon_t_z_k = reshape(cumcon_t_k, [p.nb_KFE*p.na_KFE p.nz income.ny]);
            deathin_cc_k = p.deathrate * cumcon_t_z_k(grids.loc0b0a,:,k)';
            deathin_cc_k = kron(deathin_cc_k, ones(p.nb_KFE*p.na_KFE, 1));
        end

        ind1 = 1 + p.na_KFE * p.nb_KFE * p.nz * (k-1);
        ind2 = p.na_KFE * p.nb_KFE * p.nz * k;
        RHS = reshape(c(:,:,:,k), [], 1) + ytrans_cc_k + deathin_cc_k ...
                    + cumcon_t_k(:,k) / stepsize;
        cumcon_update(ind1:ind2) = FKmats{k} * RHS;
	end
end