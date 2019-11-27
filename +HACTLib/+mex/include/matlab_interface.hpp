#ifndef _MATLAB_INTERFACE_H
#define _MATLAB_INTERFACE_H

#ifdef USING_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#endif

#include "mex.hpp"
#include "mexAdapter.hpp"


using ArrayDimensions = matlab::data::ArrayDimensions;

using Array = matlab::data::Array;

template<typename T>
using TypedArray = matlab::data::TypedArray<T>;

template<typename T>
using SparseArray = matlab::data::SparseArray<T>;

template<typename T>
using buffer_ptr_t = matlab::data::buffer_ptr_t<T>;

#include "matlab_container.hpp"

#include <vector>
#include <utility>
#include <climits>

/*
	This header file provides a MATLAB-like interface in C++ to
	MATLAB arrays passed as arguments to a mex function.
*/

typedef std::shared_ptr<matlab::engine::MATLABEngine> EnginePtr;




class MexBase : public matlab::mex::Function {
	public:

		matlab::data::ArrayFactory factory;

		EnginePtr matlabPtr = getEngine();

		void displayOnMatlab(std::ostringstream& stream);

		void displayOnMatlab(const std::string& string_output);

		void throwMatlabError(const std::string& error_msg);

		void check_same_shape(const MatlabArray<double>& arr1,
			const MatlabArray<double>& arr2);

		#if USING_EIGEN
		// Eigen::Map<Eigen::MatrixXd> to_Eigen_matrix(MatlabArray<double>& arr);

		// Eigen::Map<Eigen::VectorXd> to_Eigen_vector(MatlabArray<double> & arr);

		// Eigen::SparseMatrix<double> to_Eigen_sparse(SparseArray<double>& arr);
		#endif
};

void MexBase::displayOnMatlab(std::ostringstream& stream) {
    // Pass stream content to MATLAB fprintf function
    matlabPtr->feval(u"fprintf", 0,
        std::vector<Array>({ factory.createScalar(stream.str()) }));
    // Clear stream buffer
    stream.str("");
}

void MexBase::displayOnMatlab(const std::string& string_output) {
	std::ostringstream stream;
	stream << string_output << "\n";
	displayOnMatlab(stream);
}

void MexBase::throwMatlabError(const std::string& error_msg) {
	matlabPtr->feval(u"error", 0,
		std::vector<Array>({
			factory.createScalar(error_msg)
		}));
}

void MexBase::check_same_shape(const MatlabArray<double>& arr1,
	const MatlabArray<double>& arr2) {
	if (arr1.getDimensions() != arr2.getDimensions())
		throwMatlabError("Arrays have different shapes");
}

#if USING_EIGEN
// Eigen::Map<Eigen::MatrixXd> MexBase::to_Eigen_matrix(MatlabArray<double>& arr) {
// 	if (arr.getDimensions().size() > 2)
// 		throwMatlabError("Eigen only supports 2D arrays");

// 	return Eigen::Map<Eigen::MatrixXd>(arr.get_raw_ptr(),
// 		arr.size(0), arr.size(1));
// }

// Eigen::Map<Eigen::VectorXd> MexBase::to_Eigen_vector(MatlabArray<double>& arr) {
// 	if (arr.getDimensions().size() > 2)
// 		throwMatlabError("to_Eigen_vector requires a column vector");
// 	else if (arr.getDimensions()[1] > 1)
// 		throwMatlabError("to_Eigen_vector requires a column vector");

// 	return Eigen::Map<Eigen::VectorXd>(arr.get_raw_ptr(), arr.size(0));
// }

// Eigen::SparseMatrix<double> MexBase::to_Eigen_sparse(SparseArray<double>& arr) {
// 	std::pair<size_t, size_t> index;
// 	std::vector<double> tripletList;
// 	tripletList.reserve(arr.getNumberOfNonZeroElements());
// 	for (auto& el : arr) {
// 		index = arr.getIndex(el)
// 		tripletList.push_back(Eigen::Triplet<double>(
// 			index.first, index.second, *el));
// 	}

// 	Eigen::SparseMatrix<double> out;
// 	out.setFromTriplets(tripletList.begin(), tripletList.end());
// 	out.makeCompressed();

// 	return out;
// }

#endif


#endif