#ifndef _GSS_HEADER
#define _GSS_HEADER

#include <cmath>
#include <utility>

const double INV_GOLDEN_RATIO = 0.61803398874989479150;
const double INV_GOLDEN_RATIO_SQ = 0.38196601125;

template<typename F>
std::pair<double, double> gss(F f, const double lower, const double upper, const double tol=1e-8) {
	/*
	This function iterates over the objective function f using
	the golden section search method in the interval (a,b).

	Algorithm taken from Wikipedia.
	*/

	double a, b, c, d, diff, fc, fd;

	a = lower;
	b = upper;

	diff = b - a;
	c = a + diff * INV_GOLDEN_RATIO_SQ;
	d = a + diff * INV_GOLDEN_RATIO ;

	fc = -f(c);
	fd = -f(d);

	while (fabs(c - d) > tol) {
		if (fc < fd) {
			b = d;
			d = c;
			fd = fc;

			diff = diff * INV_GOLDEN_RATIO;
			c = a + diff * INV_GOLDEN_RATIO_SQ;
			fc = -f(c);
		}
		else {
			a = c;
			c = d;
			fc = fd;

			diff = diff * INV_GOLDEN_RATIO;
			d = a + diff * INV_GOLDEN_RATIO;
			fd = -f(d);
		}
	}

	if ( fc < fd )
		return std::make_pair(-fc, c);
	else
		return std::make_pair(-fd, d);
}

#endif