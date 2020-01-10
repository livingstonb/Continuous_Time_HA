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
	HACTLib.Checks.is_square_matrix(trans);
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

    yind = simulate(options, cumdist, discrete_cumtrans);
    final_values = values(yind(:,options.t_sim));
    statistics = create_table(final_values);
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

function yind = simulate(options, cumdist, discrete_cumtrans)
	% Performs the simulations.
	%
	% Returns
	% -------
	% yind : An array of indices into the markov process,
	%	indexed by observation and time.

	yind = zeros(options.n_sim, options.t_sim, 'uint16');
	draws = rand(options.n_sim, options.t_sim, 'single');

	[~, yind(:,1)] = max(draws(:,1) <= cumdist', [], 2);

	for t = 1:options.t_sim-1
		for k = 1:numel(cumdist)
			idx = yind(:,t) == k;
			[~, yind(idx,t+1)] = max(draws(idx,t+1)...
                <= discrete_cumtrans(k,:), [], 2);
		end
	end
end

function statistics = create_table(values)

	% Demeaned values
	values_d = values - mean(values);

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