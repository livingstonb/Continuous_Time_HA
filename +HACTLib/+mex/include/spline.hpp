#ifndef _SPLINE_HEADER
#define _SPLINE_HEADER

#include <exception>
#include <bits/stdc++.h>

template<typename T>
void spline(const T* x, const T* y, const int n, 
	const double yp1, const double ypn, T* y2);

template<typename T>
inline double splint(const T* xa, const T* ya, const T* y2a, const int n,
	const double x);

template<typename T>
void spline_eval_vec(const T* x, const T* y, const int n,
	const double yp1, const double ypn, const T* pts,
	const int npts, double* results, double* y2=nullptr) {
	/*
	This is a convenience function that both fits and
	evaluates the spline over a loop of points that share
	the same grid. The grid and values are passed as x and
	y, the grid size is passed as n, the query points are
	passed as pts, and the results variable must point to
	a preinitialized array of doubles.
	*/
	bool y2_passed;
	if (y2 == nullptr) {
		y2_passed = false;
		y2 = new double[n];
	}
	else
		y2_passed = true;

	spline(x, y, n, yp1, ypn, y2);

	for (int i=0; i<npts; ++i)
		results[i] = splint(x, y, y2, n, pts[i]);

	if (!y2_passed)
		delete[] y2;
}

template<typename T>
void spline(const T* x, const T* y, const int n, 
	const double yp1, const double ypn, T* y2) {
	/*
	Given a grid x with function values y, both of length n,
	this function approximates the second derivatives of the
	function, which are output to y2.

	Inputs yp1 and ypn are the first derivatives at the first
	point and the last point. If these are set to 1e30 or above,
	the second derivatives are assumed to be zero at these points.
	*/
	int i, k;
	double p, qn, sig, un;
	double* u = new double[n-1];

	if (yp1 > 0.99e30) {
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else {
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0])/(x[1]-x[0])-yp1);
	}

	for (int i=1; i<n-1; ++i) {
		sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
		p = sig * y2[i-1] + 2.0;
		y2[i] = (sig-1.0) / p;
		u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0 * u[i] / (x[i+1]-x[i-1]) - sig * u[i-1]) / p;
	}

	if (ypn > 0.99e30) {
		qn = 0.0;
		un = 0.0;
	}
	else {
		qn = 0.5;
		un = (3.0 / (x[n-1]-x[n-2])) * (ypn-(y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
	}

	y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);

	for (int k=n-2; k>=0; --k)
		y2[k] = y2[k] * y2[k+1] + u[k];

	delete[] u;
}

template<typename T>
inline double splint(const T* xa, const T* ya, const T* y2a, const int n,
	const double x) {
	/*
	Given a grid xa with function values ya, both of length n,
	this function outputs the interpolated function value at
	the point x. The input y2a is a vector of second derivatives
	computed in spline().
	*/
	int klo, khi, k;
	double h, b, a;

	klo = 1;
	khi = n;

	while ((khi - klo) > 1) {
		k = (khi+klo) >> 1;
		if (xa[k-1] > x)
			khi = k;
		else
			klo = k;
	}

	khi -= 1;
	klo -= 1;

	h = xa[khi] - xa[klo];
	if (h == 0)
		throw std::runtime_error("Bad xa input to splint routine");

	a = (xa[khi]-x) / h;
	b = (x-xa[klo]) / h;
	return a * ya[klo] + b * ya[khi] +
		((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h) / 6.0;
}

#endif
