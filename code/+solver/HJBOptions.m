classdef HJBOptions
	% Class used for declaring the options passed to KFESolver.

	properties (Constant)
		% Default property values.

		default_implicit (1,1) logical = false;

		default_delta (1,1) double {mustBePositive} = 1e5;

		default_HIS_maxiters (1,1) uint16 = 10;

		default_HIS_tol (1,1) double = 1e-5;

		default_HIS_start (1,1) uint16 = 2
	end

	properties (SetAccess=private)
		% Set to true for fully implicit updating.
		implicit (1,1) logical = solver.HJBOptions.default_implicit;

		% Step size.
		delta (1,1) double {mustBePositive} = solver.HJBOptions.default_delta;

		% Max number of iterations for the Howard Improvement tep.
		HIS_maxiters (1,1) uint16 = solver.HJBOptions.default_HIS_maxiters;

		% Tolerance for the Howard Improvement Step.
		HIS_tol (1,1) double = solver.HJBOptions.default_HIS_tol;

		% Number of HJB iterations before startin the Howard
		% Improvement step.
		HIS_start (1,1) uint16 = solver.HJBOptions.default_HIS_start;
	end

	methods
		function obj = HJBOptions(varargin)

			parser = inputParser;
			addParameter(parser, 'implicit', obj.default_implicit);
			addParameter(parser, 'delta', obj.default_delta);
			addParameter(parser, 'HIS_maxiters', obj.default_HIS_maxiters);
			addParameter(parser, 'HIS_tol', obj.default_HIS_tol);
			addParameter(parser, 'HIS_start', obj.default_HIS_start);

			parse(parser, varargin{:});

			for name = fieldnames(parser.Results)'
				obj.(name{:}) = parser.Results.(name{:});
			end
		end
	end

end