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
    
	% Validate input
	options = parse_inputs(varargin{:});
	HACTLib.Checks.is_square_matrix('simulate_markov', trans);
	assert(size(trans, 1) == size(values, 1),...
		"Inputs have inconsistent shapes");

	% Approximate discrete time transition matrix
	discrete_trans = options.stepsize * trans + speye(length(values));
	discrete_cumtrans = cumsum(discrete_trans, 2);

	% Compute stationary distribution
	edist = HACTLib.aux.stat_dist(trans');

	if options.normalize
		values = values ./ dot(values, edist);
	end
	cumdist = cumsum(edist);

    [Ey0, yind, quarterly_expectation] = simulate(options, cumdist, discrete_cumtrans, values);
    final_values = values(yind(:,options.t_sim));
    statistics = create_table(Ey0, final_values, quarterly_expectation);
end

function options = parse_inputs(varargin)
	parser = inputParser;
	addParameter(parser, 'n_sim', 1e4);
	addParameter(parser, 't_sim', 1e3);
	addParameter(parser, 'stepsize', 1e-3);
	addParameter(parser, 'normalize', false);
    addParameter(parser, 'initial_index', 0);
    addParameter(parser, 'r', 0.005);
    parse(parser, varargin{:});

	options = parser.Results;
end

function [Ey0, yind, quarterly_expectation] = simulate(options, cumdist, discrete_cumtrans, values)
	% Performs the simulations.
	%
	% Returns
	% -------
	% yind : An array of indices into the markov process,
	%	indexed by observation and time.

	yind = zeros(options.n_sim, options.t_sim, 'uint16');
    yvals = zeros(options.n_sim, options.t_sim, 'single');
	draws = rand(options.n_sim, options.t_sim, 'single');

    if options.initial_index <= 0
	   [~, yind(:,1)] = max(draws(:,1) <= cumdist', [], 2);
    else
        yind(:,1) = options.initial_index;
    end

	for t = 1:options.t_sim-1
		for k = 1:numel(cumdist)
			idx = yind(:,t) == k;
			[~, yind(idx,t+1)] = max(draws(idx,t+1)...
                <= discrete_cumtrans(k,:), [], 2);
		end
	end

    dt = 1 / options.t_sim;
    if options.initial_index > 0
        for t = 1:options.t_sim
            yvals(:,t) = exp(-options.r * (t-1) * dt) * dt * values(yind(:,t));
        end

        quarterly_expectation = mean(sum(yvals, 2));
    else
        quarterly_expectation = NaN;
    end

    Ey0 = mean(values(yind(:,1)));
    
end

function statistics = create_table(Ey0, values, quarterly_expectation)

	% Demeaned values
	values_d = values - mean(values);

	labels = {
        'initial mean'
        'expected quarterly inc'
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
        Ey0
        quarterly_expectation
    	mean(values)
    	min(values)
    	max(values)
        prctile(values, 0.1)
        prctile(values, 0.25)
        prctile(values, 0.5)
        prctile(values, 0.75)
        prctile(values, 0.9)
    	std(values)
    	std(log(values))
    	mean(values_d .^ 2)
    	mean(values_d .^ 3)
    	skewness(values)
    	kurtosis(values)
    };
    statistics = table(labels, moments, 'VariableNames',...
		{'Statistic', 'Value'});
end