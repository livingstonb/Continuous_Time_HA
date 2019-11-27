function statistics = simulate_markov(trans, values, varargin)
    % Simulates a continuous-time markov process and
    % produces descriptive statistics.
    %
    % Required Inputs
    % ---------------
    % trans : A square markov transition matrix, with
    %   row sums of zero, non-positive values on the
    %   main diagonal and non-negative values off the
    %   main diagonal.
    %
    % values : Values of the markov process.
    %
    % Optional Keyword Parameters
    % ---------------------------
    % n_sim : Sample size, default 1e4.
    %
    % t_sim : Number of periods to simulate,
    %   default 1e4.
    %
    % stepsize : Step size, default 1e-3.
    %
    % normalize : Boolean. If true, the process
    %   is normalized to have mean one. Default,
    %   false.
    %
    % Outputs
    % -------
    % statistics : A table containing descriptive
    %   statistics from the simulations.
    
	import HACT_Tools.Checks

	% Validate input
	options = parse_inputs(varargin{:});
	Checks.is_square_matrix(trans);
	assert(size(trans, 1) == size(values, 1),...
		"Inputs have inconsistent shapes");

	n_states = length(values);

	% Approximate discrete time transition matrix
	discrete_trans = options.stepsize * trans + speye(n_states);
	discrete_cumtrans = cumsum(discrete_trans, 2);

	% Compute stationary distribution
	edist = HACT_Tools.aux.stat_dist(trans');

	if options.normalize
		values = values ./ dot(values, edist);
	end
	cumdist = cumsum(edist);

	%% --------------------------------------------------------------------
    % Simulation
    % ---------------------------------------------------------------------
	yind = zeros(options.n_sim, options.t_sim, 'uint8');
	draws = rand(options.n_sim, options.t_sim, 'single');

	[~, yind(:,1)] = max(draws(:,1) <= cumdist', [], 2);

	for t = 1:options.t_sim-1
		for k = 1:n_states
			idx = yind(:,t) == k;
			[~, yind(idx,t+1)] = max(draws(idx,t+1)...
                <= discrete_cumtrans(k,:), [], 2);
		end
	end

	%% --------------------------------------------------------------------
    % Table
    % ---------------------------------------------------------------------
    final_values = values(yind(:,options.t_sim));
    final_values_demeaned = final_values - mean(final_values);
    labels = {
    	'mean'
    	'min'
    	'max'
        '10th pctile'
        '25th pctile'
        '50th pctile'
        '75th pctile'
        '90th pctile'
    	'std(x)'
    	'std(log x)'
    	'E[(x-xbar)^2]'
    	'E[(x-xbar)^3]'
    	'skewness'
    	'kurtosis'
    };
    moments = {
    	mean(final_values)
    	min(final_values)
    	max(final_values)
        prctile(final_values, 0.1)
        prctile(final_values, 0.25)
        prctile(final_values, 0.5)
        prctile(final_values, 0.75)
        prctile(final_values, 0.9)
    	std(final_values)
    	std(log(final_values))
    	mean(final_values_demeaned .^ 2)
    	mean(final_values_demeaned .^ 3)
    	skewness(final_values)
    	kurtosis(final_values)
    };
	statistics = table(labels, moments, 'VariableNames',...
		{'Statistic', 'Value'});
end

function options = parse_inputs(varargin)
	parser = inputParser;
	addParameter(parser, 'n_sim', 1e4);
	addParameter(parser, 't_sim', 1e3);
	addParameter(parser, 'stepsize', 1e-3);
	addParameter(parser, 'normalize', false);
    parse(parser, varargin{:});

	options = parser.Results;
end