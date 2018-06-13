// DFT library (double-precision version).
// Written by Alexander Mukhin.
// Public domain.

#include "dft.h"
#include <stdlib.h>

#define M_PI 3.14159265358979323846
#define M_SQRT2 1.41421356237309504880

// recursive complex DFT (for internal use)
static void
rdft (int N, double complex *x, double complex *y, int sy, const double complex *e) {
	int n; // counter
	double complex t0,t1,t2,t3; // temporary values
	// split-radix DIF
	if (N==1) {
		// trivial terminal case
		y[0] = x[0];
		return;
	}
	if (N==2) {
		// terminal case
		y[0] = x[0]+x[1];
	  	y[sy] = x[0]-x[1];
		return;
	}
	if (N==4) {
		// terminal case
		t0 = x[0]+x[2];
		t1 = x[1]+x[3];
		t2 = x[0]-x[2];
		t3 = I*(x[1]-x[3]);
		y[0] = t0+t1;
		y[sy] = t2-t3;
		y[2*sy] = t0-t1;
		y[3*sy] = t2+t3;
		return;
	}
	// recursion
	// prepare sub-transform inputs
	for (n=0; n<N/4; ++n) {
		t0 = x[n]+x[n+N/2];
		t1 = x[n+N/4]+x[n+3*N/4];
		t2 = x[n]-x[n+N/2];
		t3 = I*(x[n+N/4]-x[n+3*N/4]);
		x[n] = t0;
		x[n+N/4] = t1;
		x[n+N/2] = (t2-t3)*e[2*n];
		x[n+3*N/4] = (t2+t3)*e[2*n+1];
	}
	// call sub-transforms
	rdft(N/2,x,y,2*sy,e+N/2);
	rdft(N/4,x+N/2,y+sy,4*sy,e+3*N/4);
	rdft(N/4,x+3*N/4,y+3*sy,4*sy,e+3*N/4);
}

// complex DFT
void
dft (int N, double complex *x, double complex *y, const double complex *e) {
	rdft(N,x,y,1,e);
}

// fill the exponent vector for complex DFT (for internal use)
static void
fillexp_dft (int N, double complex *e) {
	int n; // counter
	while (N>4) {
		for (n=0; n<N/4; ++n) {
			*e++ = cexp(-2*M_PI*I*n/N);
			*e++ = cexp(-2*M_PI*I*3*n/N);
		}
		N /= 2;
	}
}

// allocate and fill the exponent vector for complex DFT
double complex *
mkexp_dft (int N) {
	double complex *e = (double complex*)malloc(sizeof(double complex)*N);
	fillexp_dft(N,e);
	return e;
}

// recursive unnormalized inverse complex DFT (for internal use)
static void
ridft (int N, double complex *x, double complex *y, int sy, const double complex *e) {
	int n; // counter
	double complex t0,t1,t2,t3; // temporary values
	// split-radix DIF
	if (N==1) {
		// trivial terminal case
		y[0] = x[0];
		return;
	}
	if (N==2) {
		// terminal case
		y[0] = x[0]+x[1];
	  	y[sy] = x[0]-x[1];
		return;
	}
	if (N==4) {
		// terminal case
		t0 = x[0]+x[2];
		t1 = x[1]+x[3];
		t2 = x[0]-x[2];
		t3 = I*(x[1]-x[3]);
		y[0] = t0+t1;
		y[sy] = t2+t3;
		y[2*sy] = t0-t1;
		y[3*sy] = t2-t3;
		return;
	}
	// recursion
	// prepare sub-transform inputs
	for (n=0; n<N/4; ++n) {
		t0 = x[n]+x[n+N/2];
		t1 = x[n+N/4]+x[n+3*N/4];
		t2 = x[n]-x[n+N/2];
		t3 = I*(x[n+N/4]-x[n+3*N/4]);
		x[n] = t0;
		x[n+N/4] = t1;
		x[n+N/2] = (t2+t3)*e[2*n];
		x[n+3*N/4] = (t2-t3)*e[2*n+1];
	}
	// call sub-transforms
	ridft(N/2,x,y,2*sy,e+N/2);
	ridft(N/4,x+N/2,y+sy,4*sy,e+3*N/4);
	ridft(N/4,x+3*N/4,y+3*sy,4*sy,e+3*N/4);
}

// unnormalized inverse complex DFT
void
idft (int N, double complex *x, double complex *y, const double complex *e) {
	ridft(N,x,y,1,e);
}

// fill the exponent vector for inverse complex DFT (for internal use)
static void
fillexp_idft (int N, double complex *e) {
	int n; // counter
	while (N>4) {
		for (n=0; n<N/4; ++n) {
			*e++ = cexp(2*M_PI*I*n/N);
			*e++ = cexp(2*M_PI*I*3*n/N);
		}
		N /= 2;
	}
}

// allocate and fill the exponent vector for inverse complex DFT
double complex *
mkexp_idft (int N) {
	double complex *e = (double complex*)malloc(sizeof(double complex)*N);
	fillexp_idft(N,e);
	return e;
}

// real DFT
void
realdft (int N, double *x, double *y, const double complex *e) {
	double complex *z,*w; // real vectors viewed as complex ones
	int n; // counter
	double complex a,b; // temporary values
	z = (double complex*)x;
	w = (double complex*)y;
	if (N==1) {
		// trivial case
		y[0] = x[0];
		return;
	}
	if (N==2) {
		// trivial case
		y[0] = x[0]+x[1];
	  	y[1] = x[0]-x[1];
		return;
	}
	// reduce to complex DFT of length N/2
	// do complex DFT
	dft(N/2,z,w,e+N/4);
	// recover real DFT
	w[0] = (y[0]+y[1])+I*(y[0]-y[1]);
	for (n=1; n<N/4; ++n) {
		a = (w[n]+conj(w[N/2-n]))/2;
		b = (w[n]-conj(w[N/2-n]))*e[n]/(2*I);
		w[n] = a+b;
		w[N/2-n] = conj(a-b);
	}
	w[N/4] = conj(w[N/4]);
}

// fill the exponent vector for real DFT (for internal use)
static void
fillexp_realdft (int N, double complex *e) {
	int n; // counter
	for (n=0; n<N/4; ++n)
		*e++ = cexp(-2*M_PI*I*n/N);
	fillexp_dft(N/2,e);
}

// allocate and fill the exponent vector for real DFT
double complex *
mkexp_realdft (int N) {
	double complex *e = (double complex*)malloc(sizeof(double complex)*3*N/4);
	fillexp_realdft(N,e);
	return e;
}

// unnormalized inverse real DFT
void
irealdft (int N, double *x, double *y, const double complex *e) {
	double complex *z,*w; // real vectors viewed as complex ones
	int n; // counter
	double complex a,b; // temporary values
	z = (double complex*)x;
	w = (double complex*)y;
	if (N==1) {
		// trivial case
		y[0] = x[0];
		return;
	}
	if (N==2) {
		// trivial case
		y[0] = x[0]+x[1];
	  	y[1] = x[0]-x[1];
		return;
	}
	// reduce to inverse complex DFT of length N/2
	// prepare complex DFT inputs
	z[0] = (x[0]+x[1])+I*(x[0]-x[1]);
	for (n=1; n<N/4; ++n) {
		a = z[n]+conj(z[N/2-n]);
		b = I*(z[n]-conj(z[N/2-n]))*e[n];
		z[n] = a+b;
		z[N/2-n] = conj(a-b);
	}
	z[N/4] = 2*conj(z[N/4]);
	// make inverse complex DFT
	idft(N/2,z,w,e+N/4);
}

// fill the exponent vector for inverse real DFT (for internal use)
static void
fillexp_irealdft (int N, double complex *e) {
	int n; // counter
	for (n=0; n<N/4; ++n)
		*e++ = cexp(2*M_PI*I*n/N);
	fillexp_idft(N/2,e);
}

// allocate and fill the exponent vector for inverse real DFT
double complex *
mkexp_irealdft (int N) {
	double complex *e = (double complex*)malloc(sizeof(double complex)*3*N/4);
	fillexp_irealdft(N,e);
	return e;
}

// DCT-2 of N reals
void
dct2 (int N, double *x, double *y, const double complex *e) {
	int n; // counter
	double c,s,a,b; // temporary values
	if (N==1) {
		// trivial case
		y[0] = 2*x[0];
		return;
	}
	// reduce to real DFT of length N
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n) {
		y[n] = x[2*n];
		y[N/2+n] = x[N-1-2*n];
	}
        // do real DFT
	realdft(N,y,x,e+N/2);
	// recover results
	for (n=1; n<N/2; ++n) {
		a = x[2*n];
		b = x[2*n+1];
		c = creal(e[n]);
		s = cimag(e[n]);
		y[n] = 2*(a*c-b*s);
		y[N-n] = 2*(-b*c-a*s);
	}
	// treat boundary cases
	y[N/2] = M_SQRT2*x[1];
	y[0] = 2*x[0];
}

// DST-2 of N reals
void
dst2 (int N, double *x, double *y, const double complex *e) {
	int n; // counter
	double c,s,a,b; // temporary values
	if (N==1) {
		// trivial case
		y[0] = -2*x[0];
		return;
	}
	// reduce to real DFT of length N
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n) {
		y[n] = x[2*n];
		y[N/2+n] = -x[N-1-2*n];
	}
        // do real DFT
	realdft(N,y,x,e+N/2);
	// recover results
	for (n=1; n<N/2; ++n) {
		a = x[2*n];
		b = x[2*n+1];
		c = creal(e[n]);
		s = cimag(e[n]);
		y[n] = 2*(-a*c+b*s);
		y[N-n] = 2*(b*c+a*s);
	}
	// treat boundary cases
	y[N/2] = -M_SQRT2*x[1];
	y[0] = -2*x[0];
}

// allocate and fill the exponent vector for the 2-nd type transforms
double complex *
mkexp_t2 (int N) {
	int n; // counter
	double complex *e = (double complex*)malloc(sizeof(double complex)*5*N/4);
	for (n=0; n<N/2; ++n)
		e[n] = cexp(-2*M_PI*I*n/(4*N));
	fillexp_realdft(N,e+N/2);
	return e;
}

// DCT-3 of N reals
void
dct3 (int N, double *x, double *y, const double complex *e) {
	int n; // counter
	double c,s,a,b; // temporary values
	if (N==1) {
		// trivial case
		y[0] = 2*x[0];
		return;
	}
	// reduce to inverse real DFT of length N
	// prepare sub-transform inputs
	for (n=1; n<N/2; ++n) {
		a = x[n];
		b = x[N-n];
		c = creal(e[n]);
		s = cimag(e[n]);
		y[2*n] = 2*(a*c-b*s);
		y[2*n+1] = 2*(-b*c-a*s);
	}
	y[0] = 2*x[0];
	y[1] = 2*M_SQRT2*x[N/2];
	// do inverse real DFT
	irealdft(N,y,x,e+N/2);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[2*n] = x[n];
		y[N-1-2*n] = x[N/2+n];
	}
}

// DST-3 of N reals
void
dst3 (int N, double *x, double *y, const double complex *e) {
	int n; // counter
	double c,s,a,b; // temporary values
	if (N==1) {
		// trivial case
		y[0] = -2*x[0];
		return;
	}
	// reduce to inverse real DFT of length N
	// prepare sub-transform inputs
	for (n=1; n<N/2; ++n) {
		a = x[n];
		b = x[N-n];
		c = creal(e[n]);
		s = cimag(e[n]);
		y[2*n] = 2*(-a*c+b*s);
		y[2*n+1] = 2*(b*c+a*s);
	}
	y[0] = -2*x[0];
	y[1] = -2*M_SQRT2*x[N/2];
	// do inverse real DFT
	irealdft(N,y,x,e+N/2);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[2*n] = x[n];
		y[N-1-2*n] = -x[N/2+n];
	}
}

// allocate and fill the exponent vector for the 3-rd type transforms
double complex *
mkexp_t3 (int N) {
	int n; // counter
	double complex *e = (double complex*)malloc(sizeof(double complex)*5*N/4);
	for (n=0; n<N/2; ++n)
		e[n] = cexp(-2*M_PI*I*n/(4*N));
	fillexp_irealdft(N,e+N/2);
	return e;
}

// DCT-4 of N reals
void
dct4 (int N, double *x, double *y, const double complex *e) {
	int n; // counter
	double complex *z,*w; // real vectors viewed as complex ones
	z = (double complex*)y;
	w = (double complex*)x;
	if (N==1) {
		// trivial case
		y[0] = 2*M_SQRT2*x[0];
		return;
	}
	// reduce to complex DFT of length N/2
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n)
		z[n] = 2*e[n]*(x[2*n]+I*x[N-1-2*n]);
	// do complex DFT
	dft(N/2,z,w,e+N/2);
	// recover results
	e += N;
	for (n=0; n<N/2; ++n) {
		y[2*n] = 2*creal(*e++*w[n]);
		y[2*n+1] = 2*creal(*e++*conj(w[N/2-1-n]));
	}
}

// DST-4 of N reals
void
dst4 (int N, double *x, double *y, const double complex *e) {
	int n; // counter
	double complex *z,*w; // real vectors viewed as complex ones
	z = (double complex*)y;
	w = (double complex*)x;
	if (N==1) {
		// trivial case
		y[0] = -2*M_SQRT2*x[0];
		return;
	}
	// reduce to complex DFT of length N/2
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n)
		z[n] = 2*e[n]*(x[2*n]-I*x[N-1-2*n]);
	// do complex DFT
	dft(N/2,z,w,e+N/2);
	// recover results
	e += N;
	for (n=0; n<N/2; ++n) {
		y[2*n] = 2*cimag(*e++*w[n]);
		y[2*n+1] = 2*cimag(*e++*conj(w[N/2-1-n]));
	}
}

// allocate and fill the exponent vector for the 4-th type transforms
double complex *
mkexp_t4 (int N) {
	int n; // counter
	double complex *e = (double complex*)malloc(sizeof(double complex)*2*N);
	for (n=0; n<N/2; ++n)
		e[n] = cexp(-2*M_PI*I*n/(2*N));
	fillexp_dft(N/2,e+N/2);
	for (n=0; n<N; ++n)
		e[N+n] = cexp(-2*M_PI*I*(2*n+1)/(8*N));
	return e;
}
