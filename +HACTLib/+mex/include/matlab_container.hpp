#ifndef _ARRAY_CONTAINER_H
#define _ARRAY_CONTAINER_H

#include <exception>
#include <memory>

/*
	---------- ARRAY CONTAINER CLASS -------------------------------------
	This header defines the ArrayContainer class, which serves as a
	container for the pointer underlying MATLAB Arrays.
*/

enum class MatlabType {
	structure,
	array
};

class MatlabContainer {
	public:
		MatlabType type;

		virtual ~MatlabContainer() {}
};

class MatlabStructure : public MatlabContainer {
	public:
		MatlabType type = MatlabType::structure;
};

template<typename T>
class MatlabArray : public MatlabContainer {
	buffer_ptr_t<T> dataPtr = std::unique_ptr<double [], void (*)(void*)>(nullptr, matlab::data::buffer_deleter_t());

	// T* dataPtr;

	ArrayDimensions dims;

	bool self_allocated = false;

	public:
		// MatlabArray() {}

		// MatlabArray(T* buffer, ArrayDimensions dimensions)
		// 	: dataPtr(buffer), dims(dimensions) {}

		MatlabArray(ArrayDimensions dimensions)
			: dims(dimensions), self_allocated(true) {
			matlab::data::ArrayFactory factory;
			int numel = getNumberOfElements();
			dataPtr = factory.createBuffer<T>(numel);
		}

		MatlabArray(Array array) {
			dims = array.getDimensions();
			TypedArray<double> temp = std::move(array);
			dataPtr = temp.release();
		}

		MatlabType type = MatlabType::array;

		int size(const int i) const {return dims[i];};

		ArrayDimensions getDimensions() const {return dims;}

		int getNumberOfElements() const {
			int i = 1;
			for (auto el : dims)
				i *= el;
			return i;
		}

		T* get_raw_ptr() {return dataPtr.get();}

		buffer_ptr_t<T> get_ptr() {return std::move(dataPtr);}

		T operator()(const int i) const {return dataPtr[i];}

		T operator()(const int i1, const int i2) const {
			return dataPtr[i1 + dims[0]*i2];
		}

		T operator()(const int i1, const int i2, const int i3) const {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3];
		}

		T operator()(const int i1, const int i2, const int i3, const int i4) const {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3 + dims[2] * i4];
		}

		T operator()(const int i1, const int i2, const int i3, const int i4, const int i5) const {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3 + dims[2] * i4 + dims[3] * i5];
		}

		T operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6) const {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3 + dims[2] * i4 + dims[3] * i5 + dims[4] * i6];
		}

		T& operator()(const int i) {return dataPtr[i];}

		T& operator()(const int i1, const int i2) {
			return dataPtr[i1 + dims[0]*i2];
		}

		T& operator()(const int i1, const int i2, const int i3) {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3];
		}

		T& operator()(const int i1, const int i2, const int i3, const int i4) {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3 + dims[2] * i4];
		}

		T& operator()(const int i1, const int i2, const int i3, const int i4, const int i5) {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3 + dims[2] * i4 + dims[3] * i5];
		}

		T& operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6) {
			return dataPtr[i1 + dims[0] * i2 + dims[1] * i3 + dims[2] * i4 + dims[3] * i5 + dims[4] * i6];
		}

		T operator[](const int i) const {return dataPtr[i];}

		T& operator[](const int i) {return dataPtr[i];}
};

#endif