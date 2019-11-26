classdef Asserts
	methods (Static)
		function has_shape(variable)
			% Throws if variable is not a real array
			% with the required shape.

			assert(isreal(variable),...
				"Input is not an array of real numbers")
			assert(size(variable) == shape,...
				"Input is not of the required shape")
		end

		function is_square_sparse_matrix(variable)
			% Throws if 'variable' is not a real, sparse,
			% nonempty square matrix.

			TypeAssert.is_square_matrix(variable);
			assert(issparse(variable), "Matrix is not sparse")
		end

		function is_square_matrix(variable)
			% Throws if 'variable' is not a real, nonempty,
			% square matrix.

			dim1dim2_equal = (size(object, 1) == size(object, 2);
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
	end
end