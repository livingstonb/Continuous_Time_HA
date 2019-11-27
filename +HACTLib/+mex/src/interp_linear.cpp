
#include "../include/matlab_interface.hpp"
#include "../include/linterp.hpp"

/*	
This function performs linear interpolation on gridded data.
Syntax:
	
	output = interp_linear(grid, values, points)

For a grid array, where grid(:, i1, ..., ik) represents the grid corresponding
to values(:, i1, ..., ik), this function performs interpolation over
the first dimension, looping over the remaining dimensions.
For k = 2, this function is a faster equivalent of the MATLAB code below:

	for i1 = 1:n1
		for i2 = 1:n2
			output(:, i1, i2) = interp1(grid(:, i1, i2), values(:, i1, i2), points(:, i1, i2));
		end
	end

-- grid and values must both be of dimension (n1, n2, ..., nk)
-- points must be of dimension (m, n2, ..., nk), where m need not equal n1
-- interpolation occurs along the first dimension only
-- output will be the same shape as 'points'

*/

class MexFunction : public MexBase {
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		if (inputs.size() < 3)
			throwMatlabError("Too few input arguments");
		else if (inputs.size() > 3)
			throwMatlabError("Too many input arguments");

		MatlabArray<double> grid(std::move(inputs[0]));
		MatlabArray<double> values(std::move(inputs[1]));
		MatlabArray<double> points(std::move(inputs[2]));
		MatlabArray<double> out(points.getDimensions());

		check_same_shape(grid, values);
		
		auto out_dims = out.getDimensions();
		auto grid_dims = grid.getDimensions();

		// check dimensions after first
		if (out_dims.size() != grid_dims.size())
			throwMatlabError("Grid and points have inconsistent shapes");
		
		for (int i=1; i<grid_dims.size(); ++i)
			if (out_dims[i] != grid_dims[i])
				throwMatlabError("Grid and points have inconsistent shapes");

		int outer_loops = out.getNumberOfElements() / out_dims[0];

		for (int k=0; k<outer_loops; ++k) {
			for (int i=0; i<out_dims[0]; ++i) {
				out[i+k*out_dims[0]] = linterp(
					grid, values, points[i+k*out_dims[0]],
					k*grid_dims[0], (k+1)*grid_dims[0]-1);
			}
		}

		outputs[0] = factory.createArrayFromBuffer(out.getDimensions(), out.get_ptr());
	}
};