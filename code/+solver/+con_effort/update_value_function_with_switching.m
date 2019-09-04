function Vn1 = update_value_function_with_switching(...
	p,income,A,util,Vn)

	Wn = max(Vn - p.chi,[],2);
	Wn = repmat(Wn,[1,p.nc,1]);
	z = Vn(:) - Wn(:);
	z = max(z,0);

	% LCP
	l = zeros(p.nb*p.nc*income.ny,1);
	u = inf * ones(p.nb*p.nc*income.ny,1);
	B = ((p.rho+p.deathrate)*p.delta_LCP+1)*speye(p.nb*p.nc*income.ny) ...
	     - p.delta_LCP * (A + kron(income.ytrans,speye(p.nb*p.nc)));
	q = - Vn(:) - p.delta_LCP * util(:) + B * Wn(:);
	z = aux.LCP(B,q,l,u,z,0);

	LCPerror = z' * (B*z + q);
	if abs(LCPerror) > 1e-5
	    error('large LCP error = %0.4e\n',LCPerror)
	end

	Vn1 = reshape(z(:)+Wn(:),[nb,nc,ny]);
	nrm = norm(Vn1(:)-Vn(:))
end