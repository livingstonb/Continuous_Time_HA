

clear
close all
load('/home/brian/Documents/pmf_smoothing.mat')

sm = HACTLib.computation.KernelSmoother('gaussian', true);

h = linspace(0.05, 1, numel(unique(x(:)))) .^ 3;

% x = log(0.1+x); % use h = 0.05
sm.set_for_cdf_smoothing(y, x, 0.1);
% sm.set_for_pmf_smoothing(y, x, 0.05);
sm.make_plots()

xlim([0, 0.2]);
ylim('auto')

% sorted_mat = sortrows([x(:), y(:)]);
% 
% x = sorted_mat(:,1);
% y = sorted_mat(:,2);
% 
% x_u = unique(x);
% y_u = zeros(size(x_u));
% for ix = 1:numel(x_u)
%     xval = x_u(ix);
%     
%     y_u(ix) = sum(y(x == xval));
% end
% 
% keep = (x_u >= 0) & (x_u <= 2);
% 
% sm = HACTLib.computation.KernelSmoother('gaussian');
% sm.set(x_u(keep), y_u(keep), 0.01);
% % sm.make_plots()
% sm.plot_cdfs()
% 
% xlim([0, 0.2]);
% ylim('auto')
