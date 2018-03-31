// DFT library (single-precision verson).
// Written by Alexander Mukhin.
// Public domain.

#include "dftf.h"
#include <stdlib.h>

#define M_PI 3.14159265358979323846
#define M_SQRT2 1.41421356237309504880

// recursive complex DFT (for internal use)
static void
rdftf (int N, complex float *x, complex float *y, int sy, const complex float *e) {
	int n; // counter
	complex float t0,t1,t2,t3; // temporary values
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
	rdftf(N/2,x,y,2*sy,e+N/2);
	rdftf(N/4,x+N/2,y+sy,4*sy,e+3*N/4);
	rdftf(N/4,x+3*N/4,y+3*sy,4*sy,e+3*N/4);
}

// complex DFT
void
dftf (int N, complex float *x, complex float *y, const complex float *e) {
	rdftf(N,x,y,1,e);
}

// fill the exponent vector for complex DFT (for internal use)
static void
fillexp_dftf (int N, complex float *e) {
	int n; // counter
	while (N>4) {
		for (n=0; n<N/4; ++n) {
			*e++ = cexpf(-2*M_PI*I*n/N);
			*e++ = cexpf(-2*M_PI*I*3*n/N);
		}
		N /= 2;
	}
}

// allocate and fill the exponent vector for complex DFT
complex float *
mkexp_dftf (int N) {
	complex float *e = (complex float*)malloc(sizeof(complex float)*N);
	fillexp_dftf(N,e);
	return e;
}

// recursive unnormalized inverse complex DFT (for internal use)
static void
ridftf (int N, complex float *x, complex float *y, int sy, const complex float *e) {
	int n; // counter
	complex float t0,t1,t2,t3; // temporary values
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
	ridftf(N/2,x,y,2*sy,e+N/2);
	ridftf(N/4,x+N/2,y+sy,4*sy,e+3*N/4);
	ridftf(N/4,x+3*N/4,y+3*sy,4*sy,e+3*N/4);
}

// unnormalized inverse complex DFT
void
idftf (int N, complex float *x, complex float *y, const complex float *e) {
	ridftf(N,x,y,1,e);
}

// fill the exponent vector for inverse complex DFT (for internal use)
static void
fillexp_idftf (int N, complex float *e) {
	int n; // counter
	while (N>4) {
		for (n=0; n<N/4; ++n) {
			*e++ = cexpf(2*M_PI*I*n/N);
			*e++ = cexpf(2*M_PI*I*3*n/N);
		}
		N /= 2;
	}
}

// allocate and fill the exponent vector for inverse complex DFT
complex float *
mkexp_idftf (int N) {
	complex float *e = (complex float*)malloc(sizeof(complex float)*N);
	fillexp_idftf(N,e);
	return e;
}

// real DFT
void
realdftf (int N, float *x, float *y, const complex float *e) {
	complex float *z,*w; // real vectors viewed as complex ones
	int n; // counter
	complex float a,b; // temporary values
	z = (complex float*)x;
	w = (complex float*)y;
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
	dftf(N/2,z,w,e+N/4);
	// recover real DFT
	w[0] = (y[0]+y[1])+I*(y[0]-y[1]);
	for (n=1; n<N/4; ++n) {
		a = (w[n]+conjf(w[N/2-n]))/2;
		b = (w[n]-conjf(w[N/2-n]))*e[n]/(2*I);
		w[n] = a+b;
		w[N/2-n] = conjf(a-b);
	}
	w[N/4] = conjf(w[N/4]);
}

// fill the exponent vector for real DFT (for internal use)
static void
fillexp_realdftf (int N, complex float *e) {
	int n; // counter
	for (n=0; n<N/4; ++n)
		*e++ = cexpf(-2*M_PI*I*n/N);
	fillexp_dftf(N/2,e);
}

// allocate and fill the exponent vector for real DFT
complex float *
mkexp_realdftf (int N) {
	complex float *e = (complex float*)malloc(sizeof(complex float)*3*N/4);
	fillexp_realdftf(N,e);
	return e;
}

// unnormalized inverse real DFT
void
irealdftf (int N, float *x, float *y, const complex float *e) {
	complex float *z,*w; // real vectors viewed as complex ones
	int n; // counter
	complex float a,b; // temporary values
	z = (complex float*)x;
	w = (complex float*)y;
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
		a = z[n]+conjf(z[N/2-n]);
		b = I*(z[n]-conjf(z[N/2-n]))*e[n];
		z[n] = a+b;
		z[N/2-n] = conjf(a-b);
	}
	z[N/4] = 2*conjf(z[N/4]);
	// make inverse complex DFT
	idftf(N/2,z,w,e+N/4);
}

// fill the exponent vector for inverse real DFT (for internal use)
static void
fillexp_irealdftf (int N, complex float *e) {
	int n; // counter
	for (n=0; n<N/4; ++n)
		*e++ = cexpf(2*M_PI*I*n/N);
	fillexp_idftf(N/2,e);
}

// allocate and fill the exponent vector for inverse real DFT
complex float *
mkexp_irealdftf (int N) {
	complex float *e = (complex float*)malloc(sizeof(complex float)*3*N/4);
	fillexp_irealdftf(N,e);
	return e;
}

// DCT-2 of N reals
void
dct2f (int N, float *x, float *y, const complex float *e) {
	int n; // counter
	float c,s,a,b; // temporary values
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
	realdftf(N,y,x,e+N/2);
	// recover results
	for (n=1; n<N/2; ++n) {
		a = x[2*n];
		b = x[2*n+1];
		c = crealf(e[n]);
		s = cimagf(e[n]);
		y[n] = 2*(a*c-b*s);
		y[N-n] = 2*(-b*c-a*s);
	}
	// treat boundary cases
	y[N/2] = M_SQRT2*x[1];
	y[0] = 2*x[0];
}

// DST-2 of N reals
void
dst2f (int N, float *x, float *y, const complex float *e) {
	int n; // counter
	float c,s,a,b; // temporary values
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
	realdftf(N,y,x,e+N/2);
	// recover results
	for (n=1; n<N/2; ++n) {
		a = x[2*n];
		b = x[2*n+1];
		c = crealf(e[n]);
		s = cimagf(e[n]);
		y[n] = 2*(-a*c+b*s);
		y[N-n] = 2*(b*c+a*s);
	}
	// treat boundary cases
	y[N/2] = -M_SQRT2*x[1];
	y[0] = -2*x[0];
}

// allocate and fill the exponent vector for the 2-nd type transforms
complex float *
mkexp_t2f (int N) {
	int n; // counter
	complex float *e = (complex float*)malloc(sizeof(complex float)*5*N/4);
	for (n=0; n<N/2; ++n)
		e[n] = cexpf(-2*M_PI*I*n/(4*N));
	fillexp_realdftf(N,e+N/2);
	return e;
}

// DCT-3 of N reals
void
dct3f (int N, float *x, float *y, const complex float *e) {
	int n; // counter
	float c,s,a,b; // temporary values
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
		c = crealf(e[n]);
		s = cimagf(e[n]);
		y[2*n] = 2*(a*c-b*s);
		y[2*n+1] = 2*(-b*c-a*s);
	}
	y[0] = 2*x[0];
	y[1] = 2*M_SQRT2*x[N/2];
	// do inverse real DFT
	irealdftf(N,y,x,e+N/2);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[2*n] = x[n];
		y[N-1-2*n] = x[N/2+n];
	}
}

// DST-3 of N reals
void
dst3f (int N, float *x, float *y, const complex float *e) {
	int n; // counter
	float c,s,a,b; // temporary values
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
		c = crealf(e[n]);
		s = cimagf(e[n]);
		y[2*n] = 2*(-a*c+b*s);
		y[2*n+1] = 2*(b*c+a*s);
	}
	y[0] = -2*x[0];
	y[1] = -2*M_SQRT2*x[N/2];
	// do inverse real DFT
	irealdftf(N,y,x,e+N/2);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[2*n] = x[n];
		y[N-1-2*n] = -x[N/2+n];
	}
}

// allocate and fill the exponent vector for the 3-rd type transforms
complex float *
mkexp_t3f (int N) {
	int n; // counter
	complex float *e = (complex float*)malloc(sizeof(complex float)*5*N/4);
	for (n=0; n<N/2; ++n)
		e[n] = cexpf(-2*M_PI*I*n/(4*N));
	fillexp_irealdftf(N,e+N/2);
	return e;
}

// DCT-4 of N reals
void
dct4f (int N, float *x, float *y, const complex float *e) {
	int n; // counter
	complex float *z,*w; // real vectors viewed as complex ones
	z = (complex float*)y;
	w = (complex float*)x;
	if (N==1) {
		// trivial case
		y[0] = 2*M_SQRT2*x[0];
		return;
	}
	// reduce to complex DFT of length N/2
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n)
		z[n] = 2*e[n]*(x[2*n]+I*x[N-1-2*n]);
	// do complex float DFT
	dftf(N/2,z,w,e+N/2);
	// recover results
	e += N;
	for (n=0; n<N/2; ++n) {
		y[2*n] = 2*crealf(*e++*w[n]);
		y[2*n+1] = 2*crealf(*e++*conjf(w[N/2-1-n]));
	}
}

// DST-4 of N reals
void
dst4f (int N, float *x, float *y, const complex float *e) {
	int n; // counter
	complex float *z,*w; // real vectors viewed as complex ones
	z = (complex float*)y;
	w = (complex float*)x;
	if (N==1) {
		// trivial case
		y[0] = -2*M_SQRT2*x[0];
		return;
	}
	// reduce to complex DFT of length N/2
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n)
		z[n] = 2*e[n]*(x[2*n]-I*x[N-1-2*n]);
	// do complex float DFT
	dftf(N/2,z,w,e+N/2);
	// recover results
	e += N;
	for (n=0; n<N/2; ++n) {
		y[2*n] = 2*cimagf(*e++*w[n]);
		y[2*n+1] = 2*cimagf(*e++*conjf(w[N/2-1-n]));
	}
}

// allocate and fill the exponent vector for the 4-th type transforms
complex float *
mkexp_t4f (int N) {
	int n; // counter
	complex float *e = (complex float*)malloc(sizeof(complex float)*2*N);
	for (n=0; n<N/2; ++n)
		e[n] = cexpf(-2*M_PI*I*n/(2*N));
	fillexp_dftf(N/2,e+N/2);
	for (n=0; n<N; ++n)
		e[N+n] = cexpf(-2*M_PI*I*(2*n+1)/(8*N));
	return e;
}
