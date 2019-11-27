#ifndef _BISECTION_HEADER
#define _BISECTION_HEADER

#include <cmath>
#include <utility>

template<typename F>
std::pair<double, double> bisection(F f, const double x1,
	const double x2, const double tol=1e-8, const int max_iters=50) {
	/*
	Bisection root-finding algorithm for a root bracketed between x1 and x2.
	Algorithm taken from Numerical Recipes in C.

	Return value is a pair consisting of the root and the function value at
	the root.
	*/

	double dx, fval, fmid, xmid, rtb;

	fval = f(x1);
	fmid = f(x2);
	if (fval * fmid >= 0)
		throw std::runtime_error("Function must take opposite signs at the endpoints");

	rtb = (fval < 0) ? (dx = x2-x1, x1) : (dx = x1-x2, x2);
	for (int j=0; j<=max_iters; ++j) {
		dx *= 0.5;
		xmid = rtb + dx;
		fmid = f(xmid);
		if (fmid <= 0)
			rtb = xmid;
		if ((fabs(dx) < tol) || (fmid == 0.0))
			return std::make_pair(rtb, fmid);
	}
	throw std::runtime_error("Max function evaluations exceeded");
}

template<typename F>
std::pair<double, double> secant(F f, const double x1,
	const double x2, const double tol=1e-8, const int max_iters=50) {
	/*
	Secant root-finding algorithm for a root bracketed between x1 and x2.
	Algorithm taken from Numerical Recipes in C.

	Return value is a pair consisting of the root and the function value at
	the root.
	*/
	double fl, fv, dx, swap, xl, rts;

	fl = f(x1);
	fv = f(x2);
	if (fabs(fl) < fabs(fv)) {
		rts = x1;
		xl = x2;
		swap = fl;
		fl = fv;
		fv = swap;
	}
	else {
		xl = x1;
		rts = x2;
	}

	for (int j=1; j<=max_iters; ++j) {
		dx = (xl-rts) * fv / (fv-fl);
		xl = rts;
		fl = fv;
		rts += dx;
		fv = f(rts);
		if ((fabs(dx) < tol) || (fv == 0.0))
			return std::make_pair(rts, fv);
	}
	throw std::runtime_error("Max function evaluations exceeded");
}


#endif