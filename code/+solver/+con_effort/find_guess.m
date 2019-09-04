function [V,c] = find_guess(p,income,grd)

	% ------------- unpack objects-----------------
	% assets
	bgrid_mat = squeeze(grd.b.matrix(:,1,:));
	nb = p.nb;

	% income
	y_mat = squeeze(income.y.matrix(:,1,:));
	ny = income.ny;
	% ----------------------------------------------

	% Initialization
	VbF = zeros(nb,ny);
    VbB = zeros(nb,ny);


    cF = NaN(nb,ny);
    sF = NaN(nb,ny);
    cB = NaN(nb,ny);
    sB = NaN(nb,ny);
    HcF = NaN(nb,ny,'single');
    c0 = NaN(nb,ny);

	c_0 = (1- p.wagetax) * y_mat + (p.r_b + p.deathrate*p.perfectannuities) * bgrid_mat + p.transfer;
    c_0 = max(c_0,y_mat(1,:));
	V_0	= 1/(p.rho + p.deathrate) ...
        * (aux.u_fn(c_0,p.riskaver) - aux.con_effort.penalty(bgrid_mat,p.penalty1,p.penalty2));
	
	Vbmin = 1e-8;
	Vn = V_0;
	for nn	= 1:p.maxit_HJB
		if nn > 1
			Vn = Vn1;
		end
		
		% derivatives wrt assets
		VbF(1:nb-1,:) = (Vn(2:nb,:)-Vn(1:nb-1,:)) ./ squeeze(grd.b.dF(1:nb-1,1,:));
	    VbF(1:nb-1,:) = max(VbF(1:nb-1,:),Vbmin);
	    VbB(2:nb,:) = (Vn(2:nb,:)-Vn(1:nb-1,:)) ./ squeeze(grd.b.dB(2:nb,1,:));
	    VbB(2:nb,:) = max(VbB(2:nb,:),Vbmin);

	    % consumption and savings from forward-differenced V
	    cF(1:nb-1,:) 	= VbF(1:nb-1,:).^(-1/p.riskaver);
	    cF(nb,:) 		= 1e-8;
	    sF(1:nb-1,:) 	= (1-p.wagetax) .* y_mat(1:nb-1,:)...
	    					+ bgrid_mat(1:nb-1,:) .* (p.r_b + p.deathrate*p.perfectannuities) + p.transfer - cF(1:nb-1,:);
	    sF(nb,:) 		= zeros(1,ny); % impose a state constraint at the top to improve stability
	    HcF(1:nb-1,:)	= aux.u_fn(cF(1:nb-1,:),p.riskaver) + VbF(1:nb-1,:) .* sF(1:nb-1,:);
	    HcF(nb,:) 		= -1e12 * ones(1,ny);

	    % consumption and savings from backward-differenced V
	    cB(2:nb,:)	= VbB(2:nb,:).^(-1/p.riskaver);
        cB(1,:)     = 1e-8;
        
	    sB(2:nb,:) 	= (1-p.wagetax) .* y_mat(2:nb,:)...
                        + bgrid_mat(2:nb,:) .* (p.r_b+ p.deathrate*p.perfectannuities)...
                        + p.transfer - cB(2:nb,:);
	    sB(1,:) = zeros(1,ny);
	    HcB = aux.u_fn(cB,p.riskaver) + VbB .* sB;

	    % no drift
	    c0 = (1-p.wagetax) .* y_mat...
            + bgrid_mat .* (p.r_b+ p.deathrate*p.perfectannuities) + p.transfer;
        c0 = max(c0,1e-8);
	    s0 = zeros(nb,ny);
	    Hc0 = aux.u_fn(c0,p.riskaver);

	     % Upwinding direction: consumption
	    IcF 			= (sF > 0) & ((sB>=0) | (HcF>=HcB)) & (HcF>=Hc0);
	    IcB 			= (sB < 0) & ((sF<=0) | (HcB>=HcF)) & (HcB>=Hc0);
	    assert(~any(IcF(:) & IcB(:)),'values in both IcF and IcB are true')
	    Ic0 			= ~(IcF | IcB); % 1 - IcF - IcB;
	    c 				= IcF .* cF + IcB .* cB + Ic0 .* c0;
	    s 				= IcF .* sF + IcB .* sB + Ic0 .* s0;
	    u				= aux.u_fn(c,p.riskaver) - aux.con_effort.penalty(bgrid_mat,p.penalty1,p.penalty2);
        
        % compute transition matrix
        xx = - min(s,0) ./ squeeze(grd.b.dB(:,1,:));
        xx(1,:) = 0;
        xx = xx(2:end);
        xx = xx(:);
        xx(end+1) = 0;

        zz = max(s,0) ./ squeeze(grd.b.dF(:,1,:));
        zz(nb,:) = 0;
        zz = zz(:);
        zz = [0;zz(1:end-1)];

        % A is built without income and death transitions
        yy = - max(s,0)./squeeze(grd.b.dF(:,1,:)) + min(s,0)./squeeze(grd.b.dB(:,1,:));
        yy = yy(:);

        A = spdiags([xx yy zz],[-1 0 1],nb*ny,nb*ny);
        %A = spdiags(yy(:),0,na*ny,na*ny) + spdiags(xx(:),-1,na*ny,na*ny) + spdiags(zz(:),1,na*ny,na*ny);

        % update value function
        Vn1 = NaN(nb,ny);
        for iy = 1:ny
        	% iy-iy block of A
        	Ai = A(1+(iy-1)*nb:iy*nb,1+(iy-1)*nb:iy*nb) + income.ytrans(iy,iy) * speye(nb);

        	% iy-iy block of B
        	Bi = (1 + p.delta_HJBguess * (p.rho + p.deathrate)) * speye(nb)...
                - p.delta_HJBguess * Ai;

        	uiy_stacked = u(:,iy);
        	Viy_stacked = Vn(:,iy);
        	indx_iy = ~ismember(1:ny,iy); % vector that pulls out all non-iy elements

        	% sum of trans(iy,~iy) * V(:,~iy)
        	Viyp_stacked = sum(reshape(repmat(income.ytrans(iy,indx_iy),[nb 1]),nb,ny-1)...
                .* reshape(Vn(:,indx_iy),nb,ny-1),2);

        	% DELTA bn = DELTA un + Vn
        	biy = p.delta_HJBguess * uiy_stacked + Viy_stacked + p.delta_HJBguess * Viyp_stacked;
        	Vn1iy_stacked = Bi \ biy;
        	Vn1(:,iy) = Vn1iy_stacked(:);
        end
        
	    dst = max(abs(Vn1(:)-Vn(:)));
        if mod(nn,25) == 0
            fprintf('\tHJB for guess, iteration = %i, distance = %e\n',nn,dst)
        end

	    if dst<p.crit_HJB
	    	fprintf('\tHJB for guess converged after %i iterations\n',nn);
		 	break
		elseif dst>10 && nn>1000
		 	% Not going to converge
		 	msgID = 'HJB:GuessNotConverging';
		    msg = 'HJB:GuessNotConverging';
		    HJBException = MException(msgID,msg);
		    throw(HJBException)
		end
    end
    
    V = Vn1;
end