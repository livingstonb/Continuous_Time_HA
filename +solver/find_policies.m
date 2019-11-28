function [policies, V_deriv_risky_asset_nodrift] = find_policies(p, income, grd, Vn)
    % computes policy functions on either the HJB or KFE grid

    % Parameters
    % ----------
    % p : a Params object
    %
    % income : an Income object
    %
    % grd : a Grid object
    %
    % Vn : the value function, shape (nb, na, nz, ny)
    %
    % Returns
    % -------
    % policies : a structure containing the consumption,
    %	saving, and deposits policy function
    %
    % V_deriv_risky_asset_nodrift : approximation of the
    %	first derivative of V for the case with no drift
    %	in the risky asset

    import HACTLib.computation.fd_firstorder

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

    if p.endogenous_labor == 1
        labordisutil = @(h) p.labor_disutility * (h .^ (1 + 1/p.frisch)) ./ (1 + 1/p.frisch);
        labordisutil1 = @(h) p.labor_disutility * (h .^ (1./p.frisch));
        labordisutil1inv = @(v) (v./p.labor_disutility) .^ p.frisch;
    end
    
    %% --------------------------------------------------------------------
	% UPWINDING FOR CONSUMPTION
	% ---------------------------------------------------------------------
	Vamin = 0;
    Vbmin = 1e-8;

	sspace_shape = [nb, na, nz, ny];

    % Derivatives illiquid assets
    Va = fd_firstorder(Vn, grd.a.dB, grd.a.dF, 2);
    Va.B(:,2:na,:,:) = max(Va.B(:,2:na,:,:), Vamin);
    Va.F(:,1:na-1,:,:) = max(Va.F(:,1:na-1,:,:), Vamin);

    % Derivatives liquid assets
    Vb = fd_firstorder(Vn, grd.b.dB, grd.b.dF, 1);
    Vb.B(2:nb,:,:,:) = max(Vb.B(2:nb,:,:,:), Vbmin);
    Vb.F(1:nb-1,:,:,:) = max(Vb.F(1:nb-1,:,:,:), Vbmin);

    % consumption and savings from forward-differenced V
    cF = utility1inv(Vb.F);
    cF(nb,:,:,:) = 0;

    % hours worked from forward-difference
    if p.endogenous_labor
        hoursF = labordisutil1inv(y_mat .* Vb.F);
        hoursF = min(hoursF, 1);
        sF = income.nety_HJB_liq(hoursF)
    else
        sF = income.nety_HJB_liq - cF;
    end
    sF(nb,:,:,:) = 0; % impose a state constraint at the top to improve stability

    HcF = utility(cF) + Vb.F .* sF;
    HcF(nb,:,:,:) = -1e12;

    if p.endogenous_labor
        HcF = HcF - labordisutil(hoursF);
    end

    validcF = cF > 0;

    % consumption and savings from backward-differenced V
    cB  = utility1inv(Vb.B);
    cB(1,:,:,:) = income.nety_HJB_liq(1,:,:,:);
    
    if p.endogenous_labor
        hoursB = labordisutil1inv(y_mat .* Vb.B);
        hoursB = min(hoursB, 1);
        sB = income.nety_HJB_liq(hoursB) - cB;
    else
        sB = income.nety_HJB_liq - cB;
    end
    sB(1,:,:,:) = 0;

    HcB = utility(cB) + Vb.B .* sB;
    if p.endogenous_labor
        HcB = HcB - labordisutil(hoursB);
    end

    validcB = cB > 0;

    % no drift
    c0 = income.nety_HJB_liq;
    s0 = zeros(nb,na,nz,ny);

    Hc0 = utility(c0);

    validc0 = c0 > 0;

     % Upwinding direction: consumption
    IcF = validcF & (sF > 0) & ((sB>=0) | ((HcF>=HcB) | ~validcB)) & ((HcF>=Hc0) | ~validc0);
    IcB = validcB & (sB < 0) & ((sF<=0) | ((HcB>=HcF) | ~validcF)) & ((HcB>=Hc0) | ~validc0);
    Ic0 = validc0 & ~(IcF | IcB);
    assert(isequal(IcF+IcB+Ic0,ones(nb,na,nz,ny,'logical')),'logicals do not sum to unity')
    c = IcF .* cF + IcB .* cB + Ic0 .* c0;
    s = IcF .* sF + IcB .* sB + Ic0 .* s0;

    u = utility(c);

    %% --------------------------------------------------------------------
	% UPWINDING FOR DEPOSITS
	% ---------------------------------------------------------------------
	adjcost = @(x) aux.AdjustmentCost.cost(x, grd.a.matrix, p);
	opt_d = @(x, y) aux.opt_deposits(x, y, grd.a.matrix, p);

    % Deposit decision
    dFB = opt_d(Va.F, Vb.B);
    dFB(:,na,:,:) = 0;
    dFB(1,1:na-1,:,:) = 0;
    HdFB = Va.F .* dFB - Vb.B .* (dFB + adjcost(dFB));
    HdFB(:,na,:,:) = -1.0e12;
    HdFB(1,1:na-1,:,:) = -1.0e12;
    validFB = (dFB > 0) & (HdFB > 0);

    dBF = opt_d(Va.B, Vb.F);
    dBF(:,1,:,:) = 0;
    dBF(nb,2:na,:,:) = 0;
    HdBF = Va.B .* dBF - Vb.F .* (dBF + adjcost(dBF));
    HdBF(:,1,:,:) = -1.0e12;
    HdBF(nb,2:na,:,:) = -1.0e12;
    validBF = (dBF <= - adjcost(dBF)) & (HdBF > 0);

    % Vb.B(1,2:na,:,:) = utility(cB(1,2:na,:,:));
    Vb.B(1,2:na,:,:) = utility1(cB(1,2:na,:,:));

    dBB = opt_d(Va.B, Vb.B);
    dBB(:,1,:,:) = 0;
    HdBB = Va.B .* dBB - Vb.B .* (dBB + adjcost(dBB));
    HdBB(:,1,:,:) = -1.0e12;
    validBB = (dBB > - adjcost(dBB)) & (dBB <= 0) & (HdBB > 0);
    
    if (p.OneAsset == 0) && (p.DealWithSpecialCase == 1)
        if (p.SDU == 1)
            error("Special case not coded for SDU")
        end
        
    	[H_special,c_special,d_special] = aux.deal_with_special_case(p,income,grd,r_b_mat,Va.B);

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

    d 	= IcFB .* dFB + IcBF .* dBF + IcBB .* dBB...
            + Ic_special .* d_special;

    %% --------------------------------------------------------------------
	% STORE POLICY FUNCTIONS/OTHER VARIABLES
	% ---------------------------------------------------------------------
    policies.c = c;
    policies.s = s;
    policies.d = d;
    policies.u = u;
    policies.bmin_consume_withdrawals = Ic_special;
    policies.bdot = s - adjcost(d);
    policies.adot = (p.r_a + p.deathrate*p.perfectannuities) * grd.a.matrix...
        + p.directdeposit .* y_mat + d;

    %% --------------------------------------------------------------------
    % FIRST DIFF OF VALUE FUNCTION FOR SDU WITH RETURNS RISK
    % ---------------------------------------------------------------------
    if (p.sigma_r > 0) && (p.OneAsset == 1)
        V_deriv_risky_asset_nodrift = utility1(c);
    elseif (p.sigma_r > 0) && (p.OneAsset == 0)
        V_deriv_risky_asset_nodrift = utility1(c) .* (1 + aux.AdjustmentCost.derivative(d, grd.a.matrix, p));
    else
        V_deriv_risky_asset_nodrift = [];
    end
