
a_lb = 0.01;
kappa0 = 0;
kappa1 = 2;
kappa2 = 1;

adjcost = HACTLib.aux.AdjustmentCost(a_lb, kappa0, kappa1, kappa2);

d = 0.01;
a = 0.5;
adjcost.compute_cost(d, a) / d
