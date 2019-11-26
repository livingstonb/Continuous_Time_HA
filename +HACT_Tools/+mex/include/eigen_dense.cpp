
#define USING_EIGEN true
#include "../include/matlab_interface.hpp"

/*	
	Casts a 2D MATLAB array into an Eigen::Map.
*/

class MexFunction : public MexBase {
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		double yp1 = 1.0e30;
		double ypn = 1.0e30;

		if (inputs.size() < 3)
			throwMatlabError("Too few input arguments");
		else if (inputs.size() > 3)
			throwMatlabError("Too many input arguments");

		MatlabArray<double> matlab_data(std::move(inputs[0]));
		MatlabArray<double> out(matlab_data.getDimensions());

		Eigen::Map<Eigen::MatrixXd> data(matlab_data.get_raw_ptr(),
			matlab_data.size(0), matlab_data.size(1));

		
		outputs[0] = factory.createArrayFromBuffer(out.getDimensions(), out.get_ptr());
	}
};
