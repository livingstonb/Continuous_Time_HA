function R_a_net = net_illiquid_returns(avalues, r_a,...
    tax_threshold, tax_midpt, make_plot)
    if nargin < 5
        make_plot = false;
    end
    
    if (tax_midpt <= tax_threshold)
        R_a_net = r_a .* avalues;
        return;
    end

%     decay_constant = find_constant(tax_midpt - tax_threshold);
    decay_constant = log(2) / (tax_midpt - tax_threshold);
    decay_fn = @(x) (1 - exp(-decay_constant .* x)) ./ decay_constant;

    a_low = min(avalues, tax_threshold);
    R_a_low = r_a .* a_low;
    
    a_high = max(avalues - tax_threshold, 0);
    R_a_high = r_a .* decay_fn(a_high);
    R_a_net = R_a_low + R_a_high;

    if make_plot
        plot(avalues, R_a_net);
        ylabel('Interest income from illiquid asset')
        xlabel('Level of illiquid asset, a')

        hold on

        min_y = min(R_a_net);
        max_y = max(R_a_net);
        plot([tax_threshold, tax_threshold], [min_y, max_y]);
        plot([tax_midpt, tax_midpt], [min_y, max_y]);
        legend('Total interest income', 'Tax begins',...
            'Marginal tax rate = 50%')
        
        figure()
        ihigh = a_high > 0;
        if any(ihigh)
            marginal_ra_fn = @(x) exp(-decay_constant .* x);
            mr_a = r_a .* (~ihigh + ihigh .* marginal_ra_fn(a_high));
            plot(avalues, mr_a)
            title('Marginal return, exponential decay after threshold')
            ylabel('Marginal return on the illiquid asset')
            xlabel('Level of illiquid asset, a')
            ylim([0, r_a + 0.01])
        end
    end
end

% function decay_constant = find_constant(z)
%     optfun = @(x) exp(-x * z) - 1 + x * z / 2;
%     decay_constant = fzero(optfun, 1/z);
% end