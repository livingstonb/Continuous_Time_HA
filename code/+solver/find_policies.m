function [policies, V_deriv_risky_asset_nodrift] = find_policies(p,income,grd,Vn)
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

    if p.SDU == 0
        utility = @(x) aux.u_fn(x, p.riskaver_fulldim);
        utility1 = @(x) x .^ (-p.riskaver_fulldim);
        utility1inv = @(x) x .^ (-1./p.riskaver_fulldim);
    else
        utility = @(x) rho_mat_adj .* aux.u_fn(x, p.invies);
        utility1 = @(x) rho_mat_adj .* x .^ (-p.invies);
        utility1inv = @(x) (x ./ rho_mat_adj) .^ (-1/p.invies);
    end

    if p.endogenousLabor == 1
        labordisutil = @(h) p.labor_disutility * (h .^ (1 + 1/p.frisch)) ./ (1 + 1/p.frisch);
        labordisutil1 = @(h) p.labor_disutility * (h .^ (1./p.frisch));
        labordisutil1inv = @(v) (v./p.labor_disutility) .^ p.frisch;
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

    if p.endogenousLabor == 0
        bdrift_without_c = (1-p.directdeposit-p.wagetax) .* y_mat...
            + grd.b.matrix .* (r_b_mat + p.deathrate*p.perfectannuities) + p.transfer;
    else
        bdrift_without_c = @(h) (1-p.directdeposit-p.wagetax) .* h .* y_mat...
            + grd.b.matrix .* (r_b_mat + p.deathrate*p.perfectannuities) + p.transfer;
    end

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
    cF(1:nb-1,:,:,:) = utility1inv(VbF(1:nb-1,:,:,:));
    cF(nb,:,:,:) 		= zeros(1,na,nz,ny);

    % hours worked from forward-difference
    if p.endogenousLabor == 0
        partial_drift = bdrift_without_c;
    else
        hoursF = labordisutil1inv(y_mat .* VbF);
        hoursF = min(hoursF, 1);

        partial_drift = bdrift_without_c(hoursF);
    end

    sF(1:nb-1,:,:,:)    = partial_drift(1:nb-1,:,:,:) - cF(1:nb-1,:,:,:);
    sF(nb,:,:,:)        = zeros(1,na,nz,ny); % impose a state constraint at the top to improve stability

    HcF(1:nb-1,:,:,:) = utility(cF(1:nb-1,:,:,:)) + VbF(1:nb-1,:,:,:) .* sF(1:nb-1,:,:,:);
    HcF(nb,:,:,:) 		= -1e12 * ones(1,na,nz,ny);

    if p.endogenousLabor == 1
        HcF = HcF - labordisutil(hoursF);
    end

    validcF         	= cF > 0;

    % consumption and savings from backward-differenced V
    

    if p.endogenousLabor == 0
        partial_drift = bdrift_without_c;
    else
        partial_drift = bdrift_without_c(hoursB);
        hoursB = labordisutil1inv(y_mat .* VbB);
        hoursB = min(hoursB, 1);
    end

    cB(2:nb,:,:,:)  = utility1inv(VbB(2:nb,:,:,:));
    cB(1,:,:,:) = partial_drift(1,:,:,:);
    
    sB(2:nb,:,:,:) = partial_drift(2:nb,:,:,:) - cB(2:nb,:,:,:);
    sB(1,:,:,:) = zeros(1,na,nz,ny);

    HcB = utility(cB) + VbB .* sB;

    if p.endogenousLabor == 1
        HcB = HcB - labordisutil(hoursB);
    end

    validcB = cB > 0;

    % no drift
    c0 = bdrift_without_c;
    s0 = zeros(nb,na,nz,ny);

    Hc0 = utility(c0);

    validc0         = c0 > 0;

     % Upwinding direction: consumption
    IcF 			= validcF & (sF > 0) & ((sB>=0) | ((HcF>=HcB) | ~validcB)) & ((HcF>=Hc0) | ~validc0);
    IcB 			= validcB & (sB < 0) & ((sF<=0) | ((HcB>=HcF) | ~validcF)) & ((HcB>=Hc0) | ~validc0);
    Ic0 			= validc0 & ~(IcF | IcB);
    assert(isequal(IcF+IcB+Ic0,ones(nb,na,nz,ny,'logical')),'logicals do not sum to unity')
    c 				= IcF .* cF + IcB .* cB + Ic0 .* c0;
    s 				= IcF .* sF + IcB .* sB + Ic0 .* s0;

    u = utility(c);

    %% --------------------------------------------------------------------
	% UPWINDING FOR DEPOSITS
	% ---------------------------------------------------------------------
    % Deposit decision
    dFB(2:nb,1:na-1,:,:) = aux.opt_deposits(VaF(2:nb,1:na-1,:,:),VbB(2:nb,1:na-1,:,:),grd.a.matrix(2:nb,1:na-1,:,:),p);
    dFB(:,na,:,:) = zeros(nb,1,nz,ny);
    dFB(1,1:na-1,:,:) = zeros(1,na-1,nz,ny);
    HdFB(2:nb,1:na-1,:,:) = VaF(2:nb,1:na-1,:,:) .* dFB(2:nb,1:na-1,:,:) ...
                            - VbB(2:nb,1:na-1,:,:) ...
                            .* (dFB(2:nb,1:na-1,:,:) + aux.adj_cost_fn(dFB(2:nb,1:na-1,:,:),grd.a.matrix(2:nb,1:na-1,:,:),p));
    HdFB(:,na,:,:) = -1.0e12 * ones(nb,1,nz,ny);
    HdFB(1,1:na-1,:,:) = -1.0e12 * ones(1,na-1,nz,ny);
    validFB = (dFB > 0) & (HdFB > 0);

    dBF(1:nb-1,2:na,:,:) = aux.opt_deposits(VaB(1:nb-1,2:na,:,:),VbF(1:nb-1,2:na,:,:),grd.a.matrix(1:nb-1,2:na,:,:),p);
    dBF(:,1,:,:) = zeros(nb,1,nz,ny);
    dBF(nb,2:na,:,:) = zeros(1,na-1,nz,ny);
    HdBF(1:nb-1,2:na,:,:) = VaB(1:nb-1,2:na,:,:) .* dBF(1:nb-1,2:na,:,:) - VbF(1:nb-1,2:na,:,:)...
                                 .* (dBF(1:nb-1,2:na,:,:) + aux.adj_cost_fn(dBF(1:nb-1,2:na,:,:),grd.a.matrix(1:nb-1,2:na,:,:),p));
    HdBF(:,1,:,:) = -1.0e12 * ones(nb,1,nz,ny);
    HdBF(nb,2:na,:,:) = -1.0e12 * ones(1,na-1,nz,ny);
    validBF = (dBF <= - aux.adj_cost_fn(dBF,grd.a.matrix,p)) & (HdBF > 0);

    VbB(1,2:na,:,:) = utility(cB(1,2:na,:,:));

    dBB(:,2:na,:,:) = aux.opt_deposits(VaB(:,2:na,:,:),VbB(:,2:na,:,:),grd.a.matrix(:,2:na,:,:),p);
    dBB(:,1,:,:) = zeros(nb,1,nz,ny);
    HdBB(:,2:na,:,:) = VaB(:,2:na,:,:) .* dBB(:,2:na,:,:)...
                        - VbB(:,2:na,:,:) .* (dBB(:,2:na,:,:)...
                        + aux.adj_cost_fn(dBB(:,2:na,:,:),grd.a.matrix(:,2:na,:,:),p));
    HdBB(:,1,:,:) = -1.0e12 * ones(nb,1,nz,ny);
    validBB = (dBB > - aux.adj_cost_fn(dBB,grd.a.matrix,p)) & (dBB <= 0) & (HdBB > 0);
    
    if (p.OneAsset == 0) && (p.DealWithSpecialCase == 1)
        if (p.SDU == 1)
            error("Special case not coded for SDU")
        end
        
    	[H_special,c_special,d_special] = aux.deal_with_special_case(p,income,grd,r_b_mat,VaB);

        Ic_special	= (H_special > HdFB | ~validFB) & (H_special > HdBF | ~validBF) ...
            & (H_special > HdBB | ~validBB) & (H_special > 0) & (d_special < 0)...
            & (grd.b.matrix == p.bmin) & (grd.a.matrix > 0);

        % replace c and s
        s_special = bdrift_without_c - c_special;
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
    policies.bdot = s - aux.adj_cost_fn(d,grd.a.matrix,p);
    policies.adot = (p.r_a + p.deathrate*p.perfectannuities) * grd.a.matrix...
        + p.directdeposit .* y_mat + d;

    %% --------------------------------------------------------------------
    % FIRST DIFF OF VALUE FUNCTION FOR SDU WITH RETURNS RISK
    % ---------------------------------------------------------------------
    if (p.sigma_r > 0) && (p.OneAsset == 1)
        V_deriv_risky_asset_nodrift = utility1(c);
    elseif (p.sigma_r > 0) && (p.OneAsset == 0)
        V_deriv_risky_asset_nodrift = utility1(c) .* (1 + aux.adj_cost_deriv(d, grd.a.matrix, p));
    else
        V_deriv_risky_asset_nodrift = [];
    end