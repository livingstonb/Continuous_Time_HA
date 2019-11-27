
function d = upwind_deposits(Vb, Va, adjcost, opt_d, ufn)

	na = size(Vb.B, 2);

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
    validBF = (dBF <= -adjcost(dBF)) & (HdBF > 0);

    dBB = opt_d(Va.B, Vb.B);
    dBB(:,1,:,:) = 0;
    HdBB = Va.B .* dBB - Vb.B .* (dBB + adjcost(dBB));
    HdBB(:,1,:,:) = -1.0e12;
    validBB = (dBB > - adjcost(dBB)) & (dBB <= 0) & (HdBB > 0);
end