
#include "../include/matlab_interface.hpp"
#include "../include/gss.hpp"

#include <cmath>
#include <utility>

/*
	This function loops over an input array, maximizing
	a given function at each point in the array, using a golden
	section search algorithm. The objective function is user
	defined as a lambda expression passed to gss().

	To use a value from an input array in the objective function,
	one option is to assign it to a double within the loop, and
	capture it in the lambda expression like so:
		value1 = variables[0](i, j);
		value2 = variables[1](i, j);
		gss([value1, value2](double z) {...}, lower, upper);
*/

class MexFunction : public MexBase {
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {

		MatlabArray<double> array(std::move(inputs[0]));
		MatlabArray<double> out(array.getDimensions());

		double value;
		double lower = 1.0e-8; // top of bracket
		double upper = 10.0; // bottom of bracket

		// auto gss_call = [lower, upper](double x) {return gss([x](double z) {return -pow(x-z, 2);}, lower, upper);};

		for (int i=0; i<array.getNumberOfElements(); ++i) {
			value = array[i];

			auto result = gss([value](double z) {
				//  ---------------------------------
				//  BELOW, DEFINE FUNCTION TO MAXIMIZE
				// ----------------------------------
				return -pow(z-value, 2);
			}, lower, upper);
			
			// auto result = gss_call(value);

			displayOnMatlab(std::string("value = ")+std::to_string(value));
			displayOnMatlab(std::to_string(result.second));

			// result.first is the maximum
			// result.second is the maximizer
			out[i] = result.second;
		}

		outputs[0] = factory.createArrayFromBuffer(out.getDimensions(), out.get_ptr());
	}
};