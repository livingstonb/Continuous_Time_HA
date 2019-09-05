function h = find_policies(p,income,grd,Vn)
    % output is the policy function 'h'
    % to construct policy functions on HJB grid, pass grd and Vn
    % to construct policy functions on KFE grid, pass grdKFE and VnKFE


	% ------------- unpack objects-----------------
	% asset grids
	bgrid_vec = grd.b.vec;
    cgrid_vec = grd.c.vec;
    cgrid_mat = grd.c.matrix;
	nb = numel(bgrid_vec);
    nc = numel(cgrid_vec);

	% income
	ny = income.ny;
    % ----------------------------------------------

	% Derivatives wrt c
    VcF = zeros(nb,nc,ny);
    VcB = zeros(nb,nc,ny);
    VcF(:,1:nc-1,:) = (Vn(:,2:nc,:)-Vn(:,1:nc-1,:)) ./ grd.c.dF(:,1:nc-1);
    VcB(:,2:nc,:) = (Vn(:,2:nc,:)-Vn(:,1:nc-1,:)) ./ grd.c.dB(:,2:nc);

    c = cgrid_mat;

   	% upwinding for h
    if p.chi1 == 0
        [hB,hF] = upwind_linear(p,VcF,VcB,c);
    elseif p.chi1 > 0
        [hB,hF] = upwind_nonlinear(p,VcF,VcB,c);
    end

    % Hamiltonian without u(), Va*adot, and asset penalty
    if strcmp(p.hdef,"cdot/c")
        HhF = - aux.con_effort.con_adj_cost(p,hF) + VcF .* hF .* c;
        HhB = - aux.con_effort.con_adj_cost(p,hB) + VcB .* hB .* c;
    elseif strcmp(p.hdef,"cdot")
        HhF = - aux.con_effort.con_adj_cost(p,hF) + VcF .* hF;
        HhB = - aux.con_effort.con_adj_cost(p,hB) + VcB .* hB;
    end
    Hh0 = 0;
    
   	IhF = (hF > 0) & ((hB >= 0) | (HhF>=HhB)) & (HhF>=Hh0);
   	IhB = (hB < 0) & ((hF <= 0) | (HhB>=HhF)) & (HhB>=Hh0);
	Ih0 = ~(IhF | IhB);

   	h = IhF .* hF + IhB .* hB;

end

function [hB,hF] = upwind_linear(p,VcF,VcB,c)
	% finds h with backward and forward upwinding
	if strcmp(p.hdef,"cdot/c")
        hF = (VcF.*c>p.chi0) .* p.hbar + (VcF.*c<-p.chi0) .* (-p.hbar);
        hB = (VcB.*c>p.chi0) .* p.hbar +(VcB.*c<-p.chi0) .* (-p.hbar);
    elseif strcmp(p.hdef,"cdot")
        hF = (VcF>p.chi0) .* p.hbar + (VcF<-p.chi0) .* (-p.hbar);
        hB = (VcB>p.chi0) .* p.hbar +(VcB<-p.chi0) .* (-p.hbar);
    end
end

function [hB,hF] = upwind_nonlinear(p,VcF,VcB,c)
	% finds h with backward and forward upwinding
	if strcmp(p.hdef,"cdot/c")
        hFneg = -((-VcF.*c - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hFpos = ((VcF.*c - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hF = (VcF.*c<-p.chi0) .* hFneg + (VcF.*c>p.chi0) .* hFpos;

        hBneg = -((-VcB.*c - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hBpos = ((VcB.*c - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hB = (VcB.*c<-p.chi0) .* hBneg + (VcB.*c>p.chi0) .* hBpos;
    elseif strcmp(p.hdef,"cdot")
        hFneg = -((-VcF - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hFpos = ((VcF - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hF = (VcF<-p.chi0) .* hFneg + (VcF>p.chi0) .* hFpos;

        hBneg = -((-VcB - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hBpos = ((VcB - p.chi0) / (p.chi1 * p.chi2) ) .^ (1/(p.chi2-1));
        hB = (VcB<-p.chi0) .* hBneg + (VcB>p.chi0) .* hBpos;
    end
end