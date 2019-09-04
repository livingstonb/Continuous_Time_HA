function Vn_update = update_value_function_with_drift(p,income,A,util,HJB,Vn)
	Vn_update = zeros(p.nb,p.nc,income.ny);

	for kk = 1:income.ny
		indices = 1+(kk-1)*(p.nb*p.nc):kk*(p.nb*p.nc);
        Ak = A(indices,indices);
        Bk = (1 + p.delta_HJB*(p.rho + p.deathrate)...
            	- p.delta_HJB*income.ytrans(kk,kk))*speye(p.nb*p.nc)...
                - p.delta_HJB*Ak;
        uk_stacked = reshape(util(:,:,kk),p.nb*p.nc,1);
        Vk_stacked = reshape(Vn(:,:,kk),p.nb*p.nc,1);
        indx_k = ~ismember(1:income.ny,kk);
        Vkp_stacked = sum(repmat(income.ytrans(kk,indx_k),p.nb*p.nc,1)...
                    	.* reshape(Vn(:,:,indx_k),p.nb*p.nc,income.ny-1),2);
        qk = p.delta_HJB*uk_stacked + Vk_stacked + p.delta_HJB*Vkp_stacked;
        Vn1k_stacked = Bk\qk;
        Vn_update(:,:,kk) = reshape(Vn1k_stacked,p.nb,p.nc,1);
    end
end