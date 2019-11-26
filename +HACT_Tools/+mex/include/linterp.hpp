#ifndef _LINTERP_HEADER
#define _LINTERP_HEADER

#include "locate.hpp"

/*
	This header defines a function which interpolates 'point' using
	the provided variables 'grid' and 'values.' The grid and values
	must be accessible by the bracket operator. The bottom and top
	indices of the grid and values are provided by the 'lower' and
	'upper' variables. This function uses linear interpolation and
	does not extrapolate. If the query point is below (above) the
	bottom (top) of the grid, then the value associated with the
	point at the bottom (top) of the grid is returned.
*/

template<typename T>
inline double linterp(const T& grid, const T& values, const double point,
	const int lower, const int upper) {
	int i = locate(grid, point, lower, upper);

	if (i <= lower)
		return values[lower];
	else if (i == upper + 1)
		return values[upper];
	else {
		double weight = (grid[i] - point) / (grid[i] - grid[i-1]);
		return weight * values[i-1] + (1 - weight) * values[i];
	}
}

#endif