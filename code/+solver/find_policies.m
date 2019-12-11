function [policies, V_deriv_risky_asset_nodrift] = find_policies(...
    p, income, grd, Vn, hours_bc)
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
    % hours_bc : Policy function for hours at the budget
    %   constraint, defined over the income grid.
    %
    % Returns
    % -------
    % policies : a structure containing the consumption,
    %	saving, and deposits policy function
    %
    % V_deriv_risky_asset_nodrift : approximation of the
    %	first derivative of V for the case with no drift
    %	in the risky asset

    na = numel(grd.a.vec);
    nb = numel(grd.b.vec);
    nz = p.nz;
    ny = numel(income.y.vec);
    y_mat = repmat(reshape(income.y.vec,[1 1 1 ny]),[nb na nz 1]);

    % Returns grid
    r_b_mat = p.r_b .* (grd.b.matrix>=0) +  p.r_b_borr .* (grd.b.matrix<0);

    % If using stoch diff utility, multiply utility by rho
    if p.SDU
        if numel(p.rhos) > 1
            rho_mat_adj = reshape(p.rhos,[1 1 numel(p.rhos)]);
        else
            rho_mat_adj = p.rho;
        end
        rho_mat_adj = rho_mat_adj + p.deathrate;
    else
    	rho_mat_adj = 1;
    end

    import HACTLib.model_objects.Preferences
    prefs = Preferences();
    prefs.set_crra(p.invies);

    if p.endogenous_labor
    	prefs.set_frisch(p.labordisutil, p.frisch);
    	if strcmp(grd.gtype, 'HJB')
	        nety_mat = (1-p.wagetax) * income.y.matrix;
	    else
	        nety_mat = (1-p.wagetax) * income.y.matrixKFE;
	    end
        hours_fn = @(Vb) hours_u1inv_bc_adjusted(Vb, prefs,...
            nety_mat, hours_bc);
	else
		hours_fn = {};
    end
    
    %% --------------------------------------------------------------------
	% UPWINDING FOR CONSUMPTION
	% ---------------------------------------------------------------------
	Vamin = 0;
    Vbmin = 1e-8;

	sspace_shape = [nb, na, nz, ny];
    import HACTLib.computation.fd_firstorder
    
    % Derivatives illiquid assets
    Va = fd_firstorder(Vn, grd.a.dB, grd.a.dF, 2);
    Va.B(:,2:na,:,:) = max(Va.B(:,2:na,:,:), Vamin);
    Va.F(:,1:na-1,:,:) = max(Va.F(:,1:na-1,:,:), Vamin);

    % Derivatives liquid assets
    Vb = fd_firstorder(Vn, grd.b.dB, grd.b.dF, 1);
    Vb.B(2:nb,:,:,:) = max(Vb.B(2:nb,:,:,:), Vbmin);
    Vb.F(1:nb-1,:,:,:) = max(Vb.F(1:nb-1,:,:,:), Vbmin);
    
    if strcmp(grd.gtype, 'HJB')
        net_income_liq = income.nety_HJB_liq;
    else
        net_income_liq = income.nety_KFE_liq;
    end

    import HACTLib.computation.upwind_consumption

    upwindB = upwind_consumption(net_income_liq, Vb.B,...
                'B', prefs, rho_mat_adj, hours_fn{:});
    upwindF = upwind_consumption(net_income_liq, Vb.F,...
                'F', prefs, rho_mat_adj, hours_fn{:});
    HcB = upwindB.H;
    HcF = upwindF.H;

    validcF = upwindF.c > 0;
    validcB = upwindB.c > 0;

    % no drift
    c0 = net_income_liq;
    s0 = zeros(nb,na,nz,ny);

    Hc0 = rho_mat_adj .* prefs.u(c0);
    validc0 = c0 > 0;

     % Upwinding direction: consumption
    IcF = validcF & (upwindF.s > 0) & ((upwindB.s>=0) | ((HcF>=HcB) | ~validcB)) & ((HcF>=Hc0) | ~validc0);
    IcB = validcB & (upwindB.s < 0) & ((upwindF.s<=0) | ((HcB>=HcF) | ~validcF)) & ((HcB>=Hc0) | ~validc0);
    Ic0 = validc0 & ~(IcF | IcB);
    assert(isequal(IcF+IcB+Ic0,ones(nb,na,nz,ny,'logical')),'logicals do not sum to unity')
    c = IcF .* upwindF.c + IcB .* upwindB.c + Ic0 .* c0;
    s = IcF .* upwindF.s + IcB .* upwindB.s + Ic0 .* s0;

    u = rho_mat_adj .* prefs.u(c);

    %% --------------------------------------------------------------------
	% UPWINDING FOR DEPOSITS
	% ---------------------------------------------------------------------
	adjcost = @(x) aux.AdjustmentCost.cost(x, grd.a.matrix, p);
	opt_d = @(x, y) aux.opt_deposits(x, y, grd.a.matrix, p);

    import HACTLib.computation.upwind_deposits

    Vb.B(1,2:na,:,:) = rho_mat_adj .* prefs.u1(upwindB.c(1,2:na,:,:));
    [d, I_specialcase] = upwind_deposits(Vb, Va, adjcost, opt_d);

    %% --------------------------------------------------------------------
	% STORE POLICY FUNCTIONS/OTHER VARIABLES
	% ---------------------------------------------------------------------
    policies.c = c;
    policies.s = s;
    policies.d = d;
    policies.u = u;
    policies.bmin_consume_withdrawals = I_specialcase;
    policies.bdot = s - adjcost(d);
    policies.adot = (p.r_a + p.deathrate*p.perfectannuities) * grd.a.matrix...
        + p.directdeposit .* y_mat + d;

    %% --------------------------------------------------------------------
    % DERIVATIVE OF VALUE FUNCTION FOR SDU WITH RETURNS RISK
    % ---------------------------------------------------------------------
    if (p.sigma_r > 0) && (p.OneAsset == 1)
        V_deriv_risky_asset_nodrift = rho_mat_adj .* prefs.u1(c);
    elseif (p.sigma_r > 0) && (p.OneAsset == 0)
        V_deriv_risky_asset_nodrift =rho_mat_adj .* prefs.u1(c)...
        	.* (1 + aux.AdjustmentCost.derivative(d, grd.a.matrix, p));
    else
        V_deriv_risky_asset_nodrift = [];
    end
end

function hours = hours_u1inv_bc_adjusted(Vb, prefs, nety_mat, hours_bc)
    hours = prefs.hours_u1inv(nety_mat .* Vb);
    hours(1,:,:,:) = hours_bc;
end