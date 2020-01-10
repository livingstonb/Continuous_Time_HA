#ifndef _LOCATE_HEADER
#define _LOCATE_HEADER

template<typename T>
inline int locate(const T& grid, const double entry, const int lower, const int upper) {
	int i, j, mid;

	if ( entry >= grid[upper] )
		return upper+1;
	else if ( entry <= grid[lower] )
		return lower;

	i = lower - 1;
	j = upper + 1;

	while ( (j - i) > 1 ) {
		mid = (i + j) >> 1;

		if (entry >= grid[mid])
			i = mid;
		else
			j = mid;
	}

	return (entry > grid[i]) ? j : i;
}

#endif