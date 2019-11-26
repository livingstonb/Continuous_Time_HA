classdef KFEOptions
	% Class used for declaring the options passed to KFESolver.

	properties (SetAccess=private)
		% True if using iterative procedure, false
		% otherwise.
		iterative (1,1) logical = true;

		% Step size for the iterative procedure.
		delta (1,1) double {mustBePositive} = 1e3;

		% Convergence tolerance for the iterative
		% procedure.
		tol (1,1) double {mustBePositive} = 1e-8;

		% Maximum number of iterations for the
		% iterative procedure.
		maxiters (1,1) double {mustBePositive} = 1e4;

		% If true, performs a check after a number
		% of iterations. If convergence looks unlikely,
		% an error is thrown.
		intermediate_check (1,1) logical = true;
	end
	methods
		function obj = KFEOptions(varargin)

			if nargin == 1
				if isa(varargin{1}, 'Params')
					obj.setFromParams(varargin{1});
					return;
				end
			end

			parser = inputParser;
			addParameter(parser, 'delta', 1e5);
			addParameter(parser, 'tol', 1e-8);
			addParameter(parser, 'maxiters', 1e4);
			addParameter(parser, 'iterative', true);
			addParameter(parser, 'intermediate_check', true);

			parse(parser, varargin{:});

			for name = fieldnames(parser.Results)'
				obj.(name{:}) = parser.Results.(name{:});
			end
		end
	end

end