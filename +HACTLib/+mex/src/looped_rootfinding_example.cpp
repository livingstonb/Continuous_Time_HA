
#include "../include/matlab_interface.hpp"
#include "../include/root_finding.hpp"

#include <cmath>
#include <utility>

/*
	This function loops over an input array, finding the root of
	a given function at each point in the array. The objective
	function is user defined as a lambda expression passed to the
	rootfinding algorithm.

	To use a value from an input array in the objective function,
	one option is to assign it to a double within the loop, and
	capture it in the funct object like so:
		value1 = variables[0](i, j);
		value2 = variables[1](i, j);
*/

class MexFunction : public MexBase {
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		MatlabArray<double> grid(std::move(inputs[0]));
		MatlabArray<double> out(grid.getDimensions());
		double value;
		double lower = -10.0; // top of bracket
		double upper = 10.0; // bottom of bracket

		for (int i=0; i<grid.getNumberOfElements(); ++i) {
			value = grid[i];

			auto result = bisection([value](double z) {
				/* ---------------------------------
				 BELOW, DEFINE FUNCTION TO MAXIMIZE
				----------------------------------*/
				return z - value;
			}, lower, upper);

			// result.first is the root
			// result.second is the function value
			out[i] = result.first;
		}

		outputs[0] = factory.createArrayFromBuffer(out.getDimensions(), out.get_ptr());
	}
};