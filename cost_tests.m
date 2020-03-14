
clear

amin = 0;
amax = 50;



acost = HACTLib.aux.AdjustmentCost(0.25, 0.01, 1, 1);
agrid = HACTLib.model_objects.Grid.create(amin, amax,...
	0.01, 0.9, 0.2, 50);

cost_deriv = @(d, a) -0.01 - abs(d ./ a) .^ 1;

fplot(@(d) acost.compute_cost(d, amax), [0, 10])
xlabel('Deposit rate')
ylabel('Deposit cost')

% fplot(@(d) cost_deriv(d, amax), [0, 100])