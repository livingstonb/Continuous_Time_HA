

clear

a_lb = 0.01;
kappa0 = 0;
yval = 0.5;

kappa1s = 0.1:0.1:5;
kappa2s = 0.1:0.1:5;

nk1 = numel(kappa1s);
nk2 = numel(kappa2s);

cost = zeros(nk1, nk2);
d = 0.05;

for ik1 = 1:nk1
    for ik2 = 1:nk2
        kappa1 = kappa1s(ik1);
        kappa2 = kappa2s(ik2);
        adjcostobj = HACTLib.aux.AdjustmentCost(...
            a_lb, kappa0, kappa1, kappa2);
        cost(ik1,ik2) = adjcostobj.compute_cost(d, yval) / d;
    end
end

idk1 = (kappa1s == 1);
idk2 = (kappa2s == 1);
cost(idk1, idk2)