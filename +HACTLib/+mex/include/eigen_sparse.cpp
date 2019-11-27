
#define USING_EIGEN true
#include "../include/matlab_interface.hpp"

/*	
	Casts a 2D MATLAB array into an Eigen::Map.
*/

class MexFunction : public MexBase {
	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		double yp1 = 1.0e30;
		double ypn = 1.0e30;

		if (inputs.size() < 5)
			throwMatlabError("Too few input arguments");
		else if (inputs.size() > 5)
			throwMatlabError("Too many input arguments");

		// get vectors in Eigen format
		MatlabArray<double> matlab_u(std::move(inputs[1]));
		MatlabArray<double> matlab_V(std::move(inputs[2]));
		Eigen::Map<Eigen::VectorXd> u(matlab_u.get_raw_ptr(), matlab_u.size(0));
		Eigen::Map<Eigen::VectorXd> V(matlab_V.get_raw_ptr(), matlab_V.size(0));


		// get A matrix as a Sparse Eigen matrix
		SparseArray<double> matlab_A(std::move(inputs[0]));
		std::pair<size_t, size_t> index;
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(matlab_A.getNumberOfNonZeroElements());
		for (auto iter=matlab_A.begin(); iter != matlab_A.end(); ++iter) {
			index = matlab_A.getIndex(iter);
			tripletList.push_back(Eigen::Triplet<double>(
				index.first, index.second, *iter));
		}

		Eigen::SparseMatrix<double> A(u.rows(), u.rows());
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		A.makeCompressed();

		// get constants
		TypedArray<double> matlab_rho = std::move(inputs[3]);
		TypedArray<double> matlab_delta = std::move(inputs[4]);

		double rho = matlab_rho[0];
		double delta = matlab_delta[0];

		// get output in appropriate format
		MatlabArray<double> out(matlab_V.getDimensions());
		Eigen::Map<Eigen::VectorXd> Vnew(out.get_raw_ptr(), out.size(0));

		// perform computations
		Eigen::VectorXd rhs = u + (A.triangularView<Eigen::StrictlyLower>() * V);

		Eigen::SparseMatrix<double> identity(V.rows(), V.rows());
		tripletList.clear();
		tripletList.reserve(V.rows());
		for (int i=0; i<V.rows(); ++i)
			tripletList.push_back(Eigen::Triplet<double>(i, i, 1.0));

		identity.setFromTriplets(tripletList.begin(), tripletList.end());

		A = (1 + delta * rho) * identity - delta * A.triangularView<Eigen::Upper>();
		Vnew = A.triangularView<Eigen::Upper>().solve(rhs);
		
		// return values
		outputs[0] = factory.createArrayFromBuffer(out.getDimensions(), out.get_ptr());
	}
};
