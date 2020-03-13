function R_a_net = net_illiquid_returns(avalues, r_a,...
    tax_threshold, tax_midpt, make_plot)
    if nargin < 5
        make_plot = false;
    end
    
    if (tax_midpt <= tax_threshold)
        R_a_net = r_a .* avalues;
    end

    decay_constant = find_constant(tax_midpt - tax_threshold);
    decay_fn = @(x) (1 - exp(-decay_constant .* x)) ./ decay_constant;

    a_low = min(avalues, tax_threshold);
    R_a_low = r_a .* a_low;
    
    a_high = max(avalues - tax_threshold, 0);
    R_a_high = r_a .* decay_fn(a_high);
    R_a_net = R_a_low + R_a_high;

    if make_plot
        plot(avalues, R_a_net);
        ylabel('Marginal illiquid return')
        xlabel('Level of illiquid asset, a')

        hold on

        min_y = min(R_a_net);
        max_y = max(R_a_net);
        plot([tax_threshold, tax_threshold], [min_y, max_y]);
        plot([tax_midpt, tax_midpt], [min_y, max_y]);
        legend('r_a, marginal', 'Tax begins',...
            'Tax paid on illiq interest income above threshold = 50%')
    end
end

function decay_constant = find_constant(z)
    optfun = @(x) exp(-x * z) - 1 + x * z / 2;
    decay_constant = fzero(optfun, 1/z);
end