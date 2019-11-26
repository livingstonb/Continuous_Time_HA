classdef Checks
	% A collection of assertions, mostly for checking
	% numeric arrays.

	methods (Static)
		function has_shape(variable, shape)
			% Throws if variable is not a real array
			% with the required shape.

			assert(isreal(variable),...
				"Input is not an array of real numbers")
            for k = 1:length(shape)
                assert(size(variable, k) == shape(k),...
                    "Input is not of required shape")
            end
		end

		function have_same_shape(varargin)
			% Performs assertions that the inputs are arrays
			% of the same shape.

			if nargin < 2
				return
			end

			for k = 1:nargin
				assert(isnumeric(varargin{k}) || iscell(varargin{k}),...
					"One or more of the inputs is not an array")
			end

			shape = size(varargin{1});
			for k = 2:nargin
				assert(isequal(shape, size(varargin{k})),...
					"Inputs do not have the same shape")
			end
		end

		function is_square_sparse_matrix(variable, n)
			% Throws if 'variable' is not a real, sparse,
			% nonempty square matrix.

			import HACT_Tools.Checks

			Checks.is_square_matrix(variable);
			assert(issparse(variable), "Matrix is not sparse")

			if nargin == 2
				assert(size(variable, 1) == n,...
					"Size of matrix does not match declared size")
			end
		end

		function is_square_matrix(variable)
			% Throws if 'variable' is not a real, nonempty,
			% square matrix.

			dim1dim2_equal = (size(variable, 1) == size(variable, 2));
			if ~(ismatrix(variable) && dim1dim2_equal)
				error("Input array is not a square matrix");
			end

			if size(variable, 1) == 0
				error("At least one input is empty")
			elseif ~isreal(variable)
				error("Input is complex and must be real")
			end
		end

		function is_real_integer(variable)
			msg = "Input must be real";
			assert(isreal(variable), msg);

			ScalarAssert.is_integer(variable);
		end

		function is_integer(variable)
			msg = "Input must be an integer";
			assert(isscalar(variable), msg);
			assert(round(variable == variable), msg);
		end

		function has_attributes(object, props)
			% Checks 'object' to verify that it has at least
			% the required attributes.
		
			assert(iscell(props),...
				"Second arg to has_properties must be a cell array")

			for i = 1:numel(props)
				if ~isprop(object, props{i})
					msg = sprintf("Required attribute '%s' not found", props{i});
					error(msg)
				end
			end
		end
	end
end
