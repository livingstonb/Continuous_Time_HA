
#include "../include/matlab_interface.hpp"
#include "../include/spline.hpp"

/*	
This function performs spline interpolation on gridded data.
Syntax:
	
	output = interp_spline(grid, values, points)

For a grid array, where grid(:, i1, ..., ik) represents the grid corresponding
to values(:, i1, ..., ik), this function performs interpolation over
the first dimension, looping over the remaining dimensions.
For k = 2, this function is (roughly) a faster equivalent of the MATLAB code below:

	for i1 = 1:n1
		for i2 = 1:n2
			output(:, i1, i2) = interp1(grid(:, i1, i2), values(:, i1, i2), points(:, i1, i2), 'spline');
		end
	end

-- grid and values must both be of dimension (n1, n2, ..., nk)
-- points must be of dimension (m, n2, ..., nk), where m need not equal n1
-- interpolation occurs along the first dimension only
-- output will be the same shape as 'points'

*/

class MexFunction : public MexBase {
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		double yp1 = 1.0e30;
		double ypn = 1.0e30;

		if (inputs.size() < 3)
			throwMatlabError("Too few input arguments");
		else if (inputs.size() > 3)
			throwMatlabError("Too many input arguments");

		MatlabArray<double> grid(std::move(inputs[0]));
		MatlabArray<double> values(std::move(inputs[1]));
		MatlabArray<double> points(std::move(inputs[2]));
		MatlabArray<double> out(points.getDimensions());

		if (grid.size(0) < 3)
			throwMatlabError("Too few grid points for spline interpolation");

		check_same_shape(grid, values);

		auto out_dims = out.getDimensions();
		auto grid_dims = grid.getDimensions();

		// check dimensions after first
		if (out_dims.size() != grid_dims.size())
			throwMatlabError("Grid and points have inconsistent shapes");

		double* grid_ptr = grid.get_raw_ptr();
		double* values_ptr = values.get_raw_ptr();
		double* points_ptr = points.get_raw_ptr();
		double* out_ptr = out.get_raw_ptr();

		for (int i=1; i<grid_dims.size(); ++i)
			if (out_dims[i] != grid_dims[i])
				throwMatlabError("Grid and points have inconsistent shapes");

		int outer_loops = out.getNumberOfElements() / out_dims[0];

		// perform interpolation
		double* y2 = new double[grid_dims[0]];
		for (int k=0; k<outer_loops; ++k) {
			spline_eval_vec(grid_ptr+k*grid_dims[0], values_ptr+k*grid_dims[0], grid_dims[0],
				yp1, ypn, points_ptr+k*out_dims[0], out_dims[0], out_ptr+k*out_dims[0], y2);
		}
		delete[] y2;

		outputs[0] = factory.createArrayFromBuffer(out.getDimensions(), out.get_ptr());
	}
};
