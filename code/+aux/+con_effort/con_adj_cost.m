function cost = con_adj_cost(p,h)
    cost = p.chi0 * abs(h) + p.chi1 * abs(h).^(p.chi2);
end