function [policies, Vdiff_SDU] = find_policies(p,income,grd,Vn)
    % stores policy functions and utility in policies.c, policies.d, policies.s, and policies.u
    % to construct policy functions on HJB grid, pass grd and Vn
    % to construct policy functions on KFE grid, pass grdKFE and VnKFE

    na = numel(grd.a.vec);
    nb = numel(grd.b.vec);
    nz = p.nz;
    ny = numel(income.y.vec);
    y_mat = repmat(reshape(income.y.vec,[1 1 1 ny]),[nb na nz 1]);

    % Returns grid
    r_b_mat = p.r_b .* (grd.b.matrix>=0) +  p.r_b_borr .* (grd.b.matrix<0);

    % If using stoch diff utility, multiply utility by rho
    if p.SDU == 1
        if numel(p.rhos) > 1
            rho_mat_adj = reshape(p.rhos,[1 1 numel(p.rhos)]);
        else
            rho_mat_adj = p.rho;
        end
        rho_mat_adj = rho_mat_adj + p.deathrate;
    end

    %% --------------------------------------------------------------------
	% INITIALIZATION
	% ---------------------------------------------------------------------
    VaF = zeros(nb,na,nz,ny);
    VaB = zeros(nb,na,nz,ny);
    VbF = zeros(nb,na,nz,ny);
    VbB = zeros(nb,na,nz,ny);
    
    cF = NaN(nb,na,nz,ny);
    sF = NaN(nb,na,nz,ny);
    HcF = NaN(nb,na,nz,ny,'single');

    cB = NaN(nb,na,nz,ny);
    sB = NaN(nb,na,nz,ny);

    c0 = NaN(nb,na,nz,ny);
    s0 = NaN(nb,na,nz,ny);

    dFB = NaN(nb,na,nz,ny);
    HdFB = NaN(nb,na,nz,ny,'single');
    dBF = NaN(nb,na,nz,ny);
    HdBF = NaN(nb,na,nz,ny,'single');
    dBB = NaN(nb,na,nz,ny);
    HdBB = NaN(nb,na,nz,ny,'single');

    Vamin = 0;
    Vbmin = 1e-8;

    %% --------------------------------------------------------------------
	% UPWINDING FOR CONSUMPTION
	% ---------------------------------------------------------------------
    % Derivatives illiquid assets
    VaF(:,1:na-1,:,:) = (Vn(:,2:na,:,:)-Vn(:,1:na-1,:,:)) ./ grd.a.dF(:,1:na-1);
    VaF(:,1:na-1,:,:) = max(VaF(:,1:na-1,:,:),Vamin);
    VaB(:,2:na,:,:) = (Vn(:,2:na,:,:)-Vn(:,1:na-1,:,:)) ./ grd.a.dB(:,2:na);
    VaB(:,2:na,:,:) = max(VaB(:,2:na,:,:),Vamin);

    % Derivatives liquid assets
    VbF(1:nb-1,:,:,:) = (Vn(2:nb,:,:,:)-Vn(1:nb-1,:,:,:)) ./ grd.b.dF(1:nb-1,:);
    VbF(1:nb-1,:,:,:) = max(VbF(1:nb-1,:,:,:),Vbmin); % ensure cF is well defined
    VbB(2:nb,:,:,:) = (Vn(2:nb,:,:,:)-Vn(1:nb-1,:,:,:)) ./ grd.b.dB(2:nb,:);
    VbB(2:nb,:,:,:) = max(VbB(2:nb,:,:,:),Vbmin); % ensure cF is well defined

    % consumption and savings from forward-differenced V
    if p.SDU == 0
        cF(1:nb-1,:,:,:) 	= VbF(1:nb-1,:,:,:).^(-1/p.riskaver_fulldim);
    else
        cF(1:nb-1,:,:,:) 	= (VbF(1:nb-1,:,:,:) ./ rho_mat_adj).^(-1/p.invies);
    end
    cF(nb,:,:,:) 		= zeros(1,na,nz,ny);
    sF(1:nb-1,:,:,:) 	= (1-p.directdeposit-p.wagetax) .* y_mat(1:nb-1,:,:,:)...
                        + grd.b.matrix(1:nb-1,:,:,:) .* (r_b_mat(1:nb-1,:,:,:) ...
                        + p.deathrate*p.perfectannuities) + p.transfer - cF(1:nb-1,:,:,:);
    sF(nb,:,:,:) 		= zeros(1,na,nz,ny); % impose a state constraint at the top to improve stability

    if p.SDU == 0
        HcF(1:nb-1,:,:,:) = aux.u_fn(cF(1:nb-1,:,:,:), p.riskaver) + VbF(1:nb-1,:,:,:) .* sF(1:nb-1,:,:,:);
    else
        HcF(1:nb-1,:,:,:) = aux.u_fn(cF(1:nb-1,:,:,:), p.invies, rho_mat_adj) + VbF(1:nb-1,:,:,:) .* sF(1:nb-1,:,:,:);
    end

    HcF(nb,:,:,:) 		= -1e12 * ones(1,na,nz,ny);
    validcF         	= cF > 0;

    % consumption and savings from backward-differenced V
    if p.SDU == 0
        cB(2:nb,:,:,:)	= VbB(2:nb,:,:,:).^(-1/p.riskaver_fulldim);
    else
        cB(2:nb,:,:,:)	= (VbB(2:nb,:,:,:) ./ rho_mat_adj).^(-1/p.invies);
    end
    cB(1,:,:,:) 	= ((1-p.directdeposit)-p.wagetax) .* y_mat(1,:,:,:) ...
                        + grd.b.matrix(1,:,:,:) .* (r_b_mat(1,:,:,:)+ p.deathrate*p.perfectannuities) + p.transfer;
    sB(2:nb,:,:,:) 	= ((1-p.directdeposit)-p.wagetax) .* y_mat(2:nb,:,:,:)...
                        + grd.b.matrix(2:nb,:,:,:) .* (r_b_mat(2:nb,:,:,:) ...
                        + p.deathrate*p.perfectannuities) + p.transfer - cB(2:nb,:,:,:);
    sB(1,:,:,:) 	= zeros(1,na,nz,ny);

    if p.SDU == 0
        HcB = aux.u_fn(cB, p.riskaver) + VbB .* sB;
    else
        HcB = aux.u_fn(cB, p.invies, rho_mat_adj) + VbB .* sB;
    end

    validcB = cB > 0;

    % no drift
    c0 = (1-p.directdeposit-p.wagetax) .* y_mat + grd.b.matrix .* (r_b_mat+ p.deathrate*p.perfectannuities) + p.transfer;
    s0 = zeros(nb,na,nz,ny);

    if p.SDU == 0
        Hc0 = aux.u_fn(c0, p.riskaver);
    else
        Hc0 = aux.u_fn(c0, p.invies, rho_mat_adj);
    end

    validc0         = c0 > 0;

     % Upwinding direction: consumption
    IcF 			= validcF & (sF > 0) & ((sB>=0) | ((HcF>=HcB) | ~validcB)) & ((HcF>=Hc0) | ~validc0);
    IcB 			= validcB & (sB < 0) & ((sF<=0) | ((HcB>=HcF) | ~validcF)) & ((HcB>=Hc0) | ~validc0);
    Ic0 			= validc0 & ~(IcF | IcB);
    assert(isequal(IcF+IcB+Ic0,ones(nb,na,nz,ny,'logical')),'logicals do not sum to unity')
    c 				= IcF .* cF + IcB .* cB + Ic0 .* c0;
    s 				= IcF .* sF + IcB .* sB + Ic0 .* s0;

    if p.SDU == 0
        u = aux.u_fn(c, p.riskaver);
    else
        u = aux.u_fn(c, p.invies, rho_mat_adj);
    end

    %% --------------------------------------------------------------------
	% UPWINDING FOR DEPOSITS
	% ---------------------------------------------------------------------
    % Deposit decision
    dFB(2:nb,1:na-1,:,:) = aux.two_asset.opt_deposits(VaF(2:nb,1:na-1,:,:),VbB(2:nb,1:na-1,:,:),grd.a.matrix(2:nb,1:na-1,:,:),p);
    dFB(:,na,:,:) = zeros(nb,1,nz,ny);
    dFB(1,1:na-1,:,:) = zeros(1,na-1,nz,ny);
    HdFB(2:nb,1:na-1,:,:) = VaF(2:nb,1:na-1,:,:) .* dFB(2:nb,1:na-1,:,:) ...
                            - VbB(2:nb,1:na-1,:,:) ...
                            .* (dFB(2:nb,1:na-1,:,:) + aux.two_asset.adj_cost_fn(dFB(2:nb,1:na-1,:,:),grd.a.matrix(2:nb,1:na-1,:,:),p));
    HdFB(:,na,:,:) = -1.0e12 * ones(nb,1,nz,ny);
    HdFB(1,1:na-1,:,:) = -1.0e12 * ones(1,na-1,nz,ny);
    validFB = (dFB > 0) & (HdFB > 0);

    dBF(1:nb-1,2:na,:,:) = aux.two_asset.opt_deposits(VaB(1:nb-1,2:na,:,:),VbF(1:nb-1,2:na,:,:),grd.a.matrix(1:nb-1,2:na,:,:),p);
    dBF(:,1,:,:) = zeros(nb,1,nz,ny);
    dBF(nb,2:na,:,:) = zeros(1,na-1,nz,ny);
    HdBF(1:nb-1,2:na,:,:) = VaB(1:nb-1,2:na,:,:) .* dBF(1:nb-1,2:na,:,:) - VbF(1:nb-1,2:na,:,:)...
                                 .* (dBF(1:nb-1,2:na,:,:) + aux.two_asset.adj_cost_fn(dBF(1:nb-1,2:na,:,:),grd.a.matrix(1:nb-1,2:na,:,:),p));
    HdBF(:,1,:,:) = -1.0e12 * ones(nb,1,nz,ny);
    HdBF(nb,2:na,:,:) = -1.0e12 * ones(1,na-1,nz,ny);
    validBF = (dBF <= - aux.two_asset.adj_cost_fn(dBF,grd.a.matrix,p)) & (HdBF > 0);

    if p.SDU == 0
        VbB(1,2:na,:,:) = aux.u_fn(cB(1,2:na,:,:),p.riskaver);
    else
        VbB(1,2:na,:,:) = aux.u_fn(cB(1,2:na,:,:), p.invies, rho_mat_adj);
    end
    dBB(:,2:na,:,:) = aux.two_asset.opt_deposits(VaB(:,2:na,:,:),VbB(:,2:na,:,:),grd.a.matrix(:,2:na,:,:),p);
    dBB(:,1,:,:) = zeros(nb,1,nz,ny);
    HdBB(:,2:na,:,:) = VaB(:,2:na,:,:) .* dBB(:,2:na,:,:)...
                        - VbB(:,2:na,:,:) .* (dBB(:,2:na,:,:)...
                        + aux.two_asset.adj_cost_fn(dBB(:,2:na,:,:),grd.a.matrix(:,2:na,:,:),p));
    HdBB(:,1,:,:) = -1.0e12 * ones(nb,1,nz,ny);
    validBB = (dBB > - aux.two_asset.adj_cost_fn(dBB,grd.a.matrix,p)) & (dBB <= 0) & (HdBB > 0);
    
    if (p.OneAsset == 0) && (p.DealWithSpecialCase == 1)
        if (p.SDU == 1)
            error("Special case not coded for SDU")
        end
        
    	[H_special,c_special,d_special] = aux.deal_with_special_case(p,income,grd,r_b_mat,VaB);

        Ic_special	= (H_special > HdFB | ~validFB) & (H_special > HdBF | ~validBF) ...
            & (H_special > HdBB | ~validBB) & (H_special > 0) & (d_special < 0)...
            & (grd.b.matrix == p.bmin) & (grd.a.matrix > 0);

        % replace c and s
        s_special = ((1-p.directdeposit)-p.wagetax) .* y_mat...
            + grd.b.matrix .* (r_b_mat + p.deathrate*p.perfectannuities)...
            + p.transfer - c_special;
        c(Ic_special) = c_special(Ic_special);
        s(Ic_special) = s_special(Ic_special);
        u = aux.u_fn(c,p.riskaver);
    else
        Ic_special = false(nb,na,nz,ny);
        d_special = zeros(nb,na,nz,ny);
    end
    IcFB 	= validFB & (~validBF | (HdFB >= HdBF)) ...
                & (~validBB | (HdFB >= HdBB)) & (~Ic_special);
    IcBF 	= validBF & (~validFB | (HdBF >= HdFB)) ...
                & (~validBB | (HdBF >= HdBB)) & (~Ic_special);
    IcBB 	= (~validFB | (HdBB >= HdFB)) & (~validBF | (HdBB >= HdBF)) ...
                & validBB & (~Ic_special);
    Ic00 	= (~validFB) & (~validBF) & (~validBB) & (~Ic_special);

    Isum = Ic_special+IcFB+IcBF+IcBB+Ic00;
    assert(isequal(Isum,ones(nb,na,nz,ny,'logical')),'logicals do not sum to unity')

    d 	= IcFB .* dFB + IcBF .* dBF + IcBB .* dBB + Ic00 .* zeros(nb,na,nz,ny) ...
            + Ic_special .* d_special;

    %% --------------------------------------------------------------------
	% STORE POLICY FUNCTIONS/OTHER VARIABLES
	% ---------------------------------------------------------------------
    policies.c = c;
    policies.s = s;
    policies.d = d;
    policies.u = u;
    policies.bmin_consume_withdrawals = Ic_special;
    policies.bdot = s - aux.two_asset.adj_cost_fn(d,grd.a.matrix,p);
    policies.adot = (p.r_a + p.deathrate*p.perfectannuities) * grd.a.matrix...
        + p.directdeposit .* y_mat + d;

    %% --------------------------------------------------------------------
    % FIRST DIFF OF VALUE FUNCTION FOR SDU WITH RETURNS RISK
    % ---------------------------------------------------------------------
    if (p.sigma_r > 0) && (p.OneAsset == 1)
        if p.SDU == 1
            Vb0 = rho_mat_adj .* ( c .^ (-p.invies) );
        else
            Vb0 = c .^ (-p.riskaver);
        end
        Vdiff_SDU = IcB .* VbB + IcF .* VbF + Ic0 .* Vb0;
    elseif (p.sigma_r > 0) && (p.OneAsset == 0)
        if p.SDU == 1
            Va0 = rho_mat_adj .* c .^ (-p.invies) .* (1 + aux.two_asset.adj_cost_deriv(d, grd.a.matrix, p));
        else
            Va0 = c .^ (-p.invies) .* (1 + aux.two_asset.adj_cost_deriv(d, grd.a.matrix, p));
        end
        % Vdiff_SDU = (policies.adot < 0) .* VaB + (policies.adot > 0) .* VaF...
        %     + (policies.adot == 0) .* Va0;
        Vdiff_SDU = Va0;
        % Vdiff_SDU = IcFB .* VaF + (IcBF + IcBB) .* VaB + Ic00 .* Va0;
    else
        Vdiff_SDU = [];
    end
