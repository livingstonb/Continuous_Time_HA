function Vn_update = update_value_function(p,income,A,util,HJB,Vn)
	Vn_update = zeros(p.nb*p.nc*p.nz,income.ny);

	Vn_k = reshape(Vn,[],income.ny);
	u_k = reshape(util,[],income.ny);

	for kk = 1:income.ny
		indices = 1+(kk-1)*(p.nb*p.nc*p.nz):kk*(p.nb*p.nc*p.nz);
        Ak = A(indices,indices);
        Bk = (1 + p.delta_HJB*(p.rho + p.deathrate)...
            	- p.delta_HJB*income.ytrans(kk,kk))*speye(p.nb*p.nc*p.nz)...
                - p.delta_HJB*Ak;
        indx_k = ~ismember(1:income.ny,kk);
        Vkp_stacked = sum(repmat(income.ytrans(kk,indx_k),p.nb*p.nc*p.nz,1)...
                    	.* Vn_k(:,indx_k),2);
        qk = p.delta_HJB * u_k(:,kk) + Vn_k(:,kk) + p.delta_HJB * Vkp_stacked;
        Vn_update(:,kk) = Bk \ qk;
    end

    Vn_update = reshape(Vn_update,[p.nb p.nc p.nz income.ny]);
end