classdef HJBOptions < handle
	% Class used for declaring the options passed to KFESolver.

	properties (Constant)
		% Default values
		defaults = struct(...
						'implicit', false,...
						'delta', 1e5,...
						'HIS_maxiters', 0,...
						'HIS_tol', 1e-5,...
						'HIS_start', 2...
					);
	end

	properties (SetAccess=private)
		% Set to true for fully implicit updating.
		implicit (1,1) logical = obj.defaults.implicit;

		% Step size.
		delta (1,1) double {mustBePositive} = obj.defaults.delta;

		% Max number of iterations for the Howard Improvement tep.
		HIS_maxiters (1,1) uint16 = obj.defaults.HIS_maxiters;

		% Tolerance for the Howard Improvement Step.
		HIS_tol (1,1) double = obj.defaults.HIS_tol;

		% Number of HJB iterations before startin the Howard
		% Improvement step.
		HIS_start (1,1) uint16 = obj.defaults.HIS_start;
	end

	methods
		function obj = HJBOptions(varargin)
			% Creates an HJBOptions object via kew
			options = HACT_Tools.aux.parse_keyvalue_pairs(...
				obj.defaults, varargin{:});
			obj.set_from_struct(options);
		end

		function set_from_struct(options)
			for name = fieldnames(options)'
				obj.(name{1}) = options.(name{1});
			end
		end
	end

end