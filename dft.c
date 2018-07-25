// DFT library (double-precision version).
// Written by Alexander Mukhin.
// Public domain.

#include "dft.h"
#include <stdlib.h>

#define M_PI 3.14159265358979323846
#define M_SQRT2 1.41421356237309504880

// *** meta-functions ***

// make a strided any-dimensional complex transform
// by repeated application of its strided one-dimensional routine
inline static void
make_complex_transform (
	double complex *x, // source data
	double complex *y, // destination buffer
	int sy, // stride on y
	const struct aux *a, // aux data
	void (*s_1d) \
	(double complex*,double complex*,int,const struct aux*) // strided 1d xform
) {
	if (a->sub2==NULL)
		(*s_1d)(x,y,sy,a);
	else {
		int N1=a->sub1->N,N2=a->sub2->N; // transform lengths
		int n; // counter
		double complex *t=a->t; // temporary buffer
		for (n=0; n<N2; ++n)
			make_complex_transform(x+n*N1,t+n,N2,a->sub1,s_1d);
		for (n=0; n<N1; ++n)
			(*s_1d)(t+n*N2,y+sy*n,sy*N1,a->sub2);
	}
}

// make a strided any-dimensional real transform
// by repeated application of its strided one-dimensional routine
inline static void
make_real_transform (
	double *x, // source data
	double *y, // destination buffer
	int sy, // stride on y
	const struct aux *a, // aux data
	void (*s_1d) \
	(double*,double*,int,const struct aux*) // strided 1d xform
) {
	if (a->sub2==NULL)
		(*s_1d)(x,y,sy,a);
	else {
		int N1=a->sub1->N,N2=a->sub2->N; // transform lengths
		int n; // counter
		double *t=a->t; // temporary buffer
		for (n=0; n<N2; ++n)
			make_real_transform(x+n*N1,t+n,N2,a->sub1,s_1d);
		for (n=0; n<N1; ++n)
			(*s_1d)(t+n*N2,y+sy*n,sy*N1,a->sub2);
	}
}

// *** complex transforms ***

// recursive strided one-dimensional DFT
inline static void
rs_dft_1d (
	int N, // transform length
	double complex *x, // source data
	double complex *t, // temporary buffer
	double complex *y, // destination buffer
	int sy, // stride on y
	const double complex *e // exponent vector
) {
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
		t0 = x[0]+x[1];
	  	t1 = x[0]-x[1];
		y[0] = t0;
	  	y[sy] = t1;
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
	if (N==8) {
		// terminal case
		double complex t00,t01,t02,t03;
		double complex t10,t11,t12,t13;
		const double complex E1=0.70710678118654752440*(1-I);
		const double complex E3=0.70710678118654752440*(-1-I);
		t0 = x[0]+x[4];
		t1 = x[2]+x[6];
		t2 = x[0]-x[4];
		t3 = I*(x[2]-x[6]);
		t00 = t0+t1;
		t01 = t2-t3;
		t02 = t0-t1;
		t03 = t2+t3;
		t0 = x[1]+x[5];
		t1 = x[3]+x[7];
		t2 = x[1]-x[5];
		t3 = I*(x[3]-x[7]);
		t10 = t0+t1;
		t11 = (t2-t3)*E1;
		t12 = (t0-t1)*(-I);
		t13 = (t2+t3)*E3;
		y[0] = t00+t10;
		y[sy] = t01+t11;
		y[2*sy] = t02+t12;
		y[3*sy] = t03+t13;
		y[4*sy] = t00-t10;
		y[5*sy] = t01-t11;
		y[6*sy] = t02-t12;
		y[7*sy] = t03-t13;
		return;
	}
	// recursion
	// prepare sub-transform inputs
	for (n=0; n<N/4; ++n) {
		t0 = x[n]+x[n+N/2];
		t1 = x[n+N/4]+x[n+3*N/4];
		t2 = x[n]-x[n+N/2];
		t3 = I*(x[n+N/4]-x[n+3*N/4]);
		t[n] = t0;
		t[n+N/4] = t1;
		t[n+N/2] = (t2-t3)*e[2*n];
		t[n+3*N/4] = (t2+t3)*e[2*n+1];
	}
	// call sub-transforms
	rs_dft_1d(N/2,t,t,y,2*sy,e+N/2);
	rs_dft_1d(N/4,t+N/2,t+N/2,y+sy,4*sy,e+3*N/4);
	rs_dft_1d(N/4,t+3*N/4,t+3*N/4,y+3*sy,4*sy,e+3*N/4);
}

// strided one-dimensional DFT
inline static void
s_dft_1d (double complex *x, double complex *y, int sy, const struct aux *a) {
	rs_dft_1d(a->N,x,a->t,y,sy,a->e);
}

// strided DFT of arbitrary dimension
inline static void
s_dft (double complex *x, double complex *y, int sy, const struct aux *a) {
	make_complex_transform(x,y,sy,a,s_dft_1d);
}

// user interface
void
dft (double complex *x, double complex *y, const struct aux *a) {
	s_dft(x,y,1,a);
}

// recursive strided one-dimensional inverse DFT
inline static void
rs_idft_1d (
	int N, // transform length
	double complex *x, // source data
	double complex *t, // temporary buffer
	double complex *y, // destination buffer
	int sy, // stride on y
	const double complex *e // exponent vector
) {
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
		t0 = x[0]+x[1];
	  	t1 = x[0]-x[1];
		y[0] = t0;
	  	y[sy] = t1;
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
	if (N==8) {
		// terminal case
		double complex t00,t01,t02,t03;
		double complex t10,t11,t12,t13;
		const double complex E1=0.70710678118654752440*(1+I);
		const double complex E3=0.70710678118654752440*(-1+I);
		t0 = x[0]+x[4];
		t1 = x[2]+x[6];
		t2 = x[0]-x[4];
		t3 = I*(x[2]-x[6]);
		t00 = t0+t1;
		t01 = t2+t3;
		t02 = t0-t1;
		t03 = t2-t3;
		t0 = x[1]+x[5];
		t1 = x[3]+x[7];
		t2 = x[1]-x[5];
		t3 = I*(x[3]-x[7]);
		t10 = t0+t1;
		t11 = (t2+t3)*E1;
		t12 = (t0-t1)*I;
		t13 = (t2-t3)*E3;
		y[0] = t00+t10;
		y[sy] = t01+t11;
		y[2*sy] = t02+t12;
		y[3*sy] = t03+t13;
		y[4*sy] = t00-t10;
		y[5*sy] = t01-t11;
		y[6*sy] = t02-t12;
		y[7*sy] = t03-t13;
		return;
	}
	// recursion
	// prepare sub-transform inputs
	for (n=0; n<N/4; ++n) {
		t0 = x[n]+x[n+N/2];
		t1 = x[n+N/4]+x[n+3*N/4];
		t2 = x[n]-x[n+N/2];
		t3 = I*(x[n+N/4]-x[n+3*N/4]);
		t[n] = t0;
		t[n+N/4] = t1;
		t[n+N/2] = (t2+t3)*conj(e[2*n]);
		t[n+3*N/4] = (t2-t3)*conj(e[2*n+1]);
	}
	// call sub-transforms
	rs_idft_1d(N/2,t,t,y,2*sy,e+N/2);
	rs_idft_1d(N/4,t+N/2,t+N/2,y+sy,4*sy,e+3*N/4);
	rs_idft_1d(N/4,t+3*N/4,t+3*N/4,y+3*sy,4*sy,e+3*N/4);
}

// strided one-dimensional inverse DFT
inline static void
s_idft_1d (double complex *x, double complex *y, int sy, const struct aux *a) {
	rs_idft_1d(a->N,x,a->t,y,sy,a->e);
}

// strided inverse DFT of arbitrary dimension
inline static void
s_idft (double complex *x, double complex *y, int sy, const struct aux *a) {
	make_complex_transform(x,y,sy,a,s_idft_1d);
}

// user interface
void
idft (double complex *x, double complex *y, const struct aux *a) {
	s_idft(x,y,1,a);
}

// *** real transforms ***

// one-dimensional real DFT
void
realdft_1d (double *x, double *y, const struct aux *a) {
	double complex *z,*w; // real vectors viewed as complex ones
	int n; // counter
	double complex u,v; // temporary values
	int N=a->N; // transform length
	double complex *e=a->e; // exponent vector
	z = (double complex*)x;
	w = (double complex*)y;
	if (N==1) {
		// trivial case
		y[0] = x[0];
		return;
	}
	if (N==2) {
		// trivial case
		double t0,t1; // temporary values
		t0 = x[0]+x[1];
	  	t1 = x[0]-x[1];
		y[0] = t0;
		y[1] = t1;
		return;
	}
	// reduce to complex DFT of length N/2
	// do complex DFT
	dft(z,w,a->sub1);
	// recover real DFT
	w[0] = (y[0]+y[1])+I*(y[0]-y[1]);
	for (n=1; n<N/4; ++n) {
		u = (w[n]+conj(w[N/2-n]))/2;
		v = (w[n]-conj(w[N/2-n]))*e[n]/(2*I);
		w[n] = u+v;
		w[N/2-n] = conj(u-v);
	}
	w[N/4] = conj(w[N/4]);
}

// one-dimensional inverse real DFT
void
irealdft_1d (double *x, double *y, const struct aux *a) {
	double complex *z,*w; // real vectors viewed as complex ones
	int n; // counter
	double complex u,v; // temporary values
	int N=a->N; // transform length
	double complex *e=a->e; // exponent vector
	z = (double complex*)x;
	w = (double complex*)y;
	if (N==1) {
		// trivial case
		y[0] = x[0];
		return;
	}
	if (N==2) {
		// trivial case
		double t0,t1; // temporary values
		t0 = x[0]+x[1];
	  	t1 = x[0]-x[1];
		y[0] = t0;
		y[1] = t1;
		return;
	}
	// reduce to inverse complex DFT of length N/2
	// prepare complex DFT inputs
	w[0] = (x[0]+x[1])+I*(x[0]-x[1]);
	for (n=1; n<N/4; ++n) {
		u = z[n]+conj(z[N/2-n]);
		v = I*(z[n]-conj(z[N/2-n]))*conj(e[n]);
		w[n] = u+v;
		w[N/2-n] = conj(u-v);
	}
	w[N/4] = 2*conj(z[N/4]);
	// make inverse complex DFT
	idft(w,w,a->sub1);
}

// *** real symmetric transforms ***

// strided one-dimensional DCT-2
inline static void
s_dct2_1d (double *x, double *y, int sy, const struct aux *a) {
	int n; // counter
	double c,s,u,v; // temporary values
	int N=a->N; // transform length
	double *t=a->t; // temporary buffer
	double complex *e=a->e; // exponent vector
	if (N==1) {
		// trivial case
		y[0] = 2*x[0];
		return;
	}
	// reduce to real DFT of length N
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n) {
		t[n] = x[2*n];
		t[N/2+n] = x[N-1-2*n];
	}
        // do real DFT in-place
	realdft_1d(t,t,a->sub1);
	// recover results
	for (n=1; n<N/2; ++n) {
		u = t[2*n];
		v = t[2*n+1];
		c = creal(e[n]);
		s = cimag(e[n]);
		y[sy*n] = 2*(u*c-v*s);
		y[sy*(N-n)] = 2*(-v*c-u*s);
	}
	// treat boundary cases
	y[sy*N/2] = M_SQRT2*t[1];
	y[0] = 2*t[0];
}

// strided DCT-2 of arbitrary dimension
inline static void
s_dct2 (double *x, double *y, int sy, const struct aux *a) {
	make_real_transform(x,y,sy,a,s_dct2_1d);
}

// user interface
void
dct2 (double *x, double *y, const struct aux *a) {
	s_dct2(x,y,1,a);
}

// strided one-dimensional DST-2
inline static void
s_dst2_1d (double *x, double *y, int sy, const struct aux *a) {
	int n; // counter
	double c,s,u,v; // temporary values
	int N=a->N; // transform length
	double *t=a->t; // temporary buffer
	double complex *e=a->e; // exponent vector
	if (N==1) {
		// trivial case
		y[0] = -2*x[0];
		return;
	}
	// reduce to real DFT of length N
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n) {
		t[n] = x[2*n];
		t[N/2+n] = -x[N-1-2*n];
	}
        // do real DFT in-place
	realdft_1d(t,t,a->sub1);
	// recover results
	for (n=1; n<N/2; ++n) {
		u = t[2*n];
		v = t[2*n+1];
		c = creal(e[n]);
		s = cimag(e[n]);
		y[sy*n] = 2*(-u*c+v*s);
		y[sy*(N-n)] = 2*(v*c+u*s);
	}
	// treat boundary cases
	y[sy*N/2] = -M_SQRT2*t[1];
	y[0] = -2*t[0];
}

// strided DST-2 of arbitrary dimension
inline static void
s_dst2 (double *x, double *y, int sy, const struct aux *a) {
	make_real_transform(x,y,sy,a,s_dst2_1d);
}

// user interface
void
dst2 (double *x, double *y, const struct aux *a) {
	s_dst2(x,y,1,a);
}

// strided one-dimensional DCT-3
inline static void
s_dct3_1d (double *x, double *y, int sy, const struct aux *a) {
	int n; // counter
	double c,s,u,v; // temporary values
	int N=a->N; // transform length
	double *t=a->t; // temporary buffer
	double complex *e=a->e; // exponent vector
	if (N==1) {
		// trivial case
		y[0] = 2*x[0];
		return;
	}
	// reduce to inverse real DFT of length N
	// prepare sub-transform inputs
	for (n=1; n<N/2; ++n) {
		u = x[n];
		v = x[N-n];
		c = creal(e[n]);
		s = cimag(e[n]);
		t[2*n] = 2*(u*c-v*s);
		t[2*n+1] = 2*(-v*c-u*s);
	}
	t[0] = 2*x[0];
	t[1] = 2*M_SQRT2*x[N/2];
	// do inverse real DFT in-place
	irealdft_1d(t,t,a->sub1);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[sy*2*n] = t[n];
		y[sy*(N-1-2*n)] = t[N/2+n];
	}
}

// strided DCT-3 of arbitrary dimension
inline static void
s_dct3 (double *x, double *y, int sy, const struct aux *a) {
	make_real_transform(x,y,sy,a,s_dct3_1d);
}

// user interface
void
dct3 (double *x, double *y, const struct aux *a) {
	s_dct3(x,y,1,a);
}

// strided one-dimensional DST-3
inline static void
s_dst3_1d (double *x, double *y, int sy, const struct aux *a) {
	int n; // counter
	double c,s,u,v; // temporary values
	int N=a->N; // transform length
	double *t=a->t; // temporary buffer
	double complex *e=a->e; // exponent vector
	if (N==1) {
		// trivial case
		y[0] = -2*x[0];
		return;
	}
	// reduce to inverse real DFT of length N
	// prepare sub-transform inputs
	for (n=1; n<N/2; ++n) {
		u = x[n];
		v = x[N-n];
		c = creal(e[n]);
		s = cimag(e[n]);
		t[2*n] = 2*(-u*c+v*s);
		t[2*n+1] = 2*(v*c+u*s);
	}
	t[0] = -2*x[0];
	t[1] = -2*M_SQRT2*x[N/2];
	// do inverse real DFT in-place
	irealdft_1d(t,t,a->sub1);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[sy*2*n] = t[n];
		y[sy*(N-1-2*n)] = -t[N/2+n];
	}
}

// strided DST-3 of arbitrary dimension
inline static void
s_dst3 (double *x, double *y, int sy, const struct aux *a) {
	make_real_transform(x,y,sy,a,s_dst3_1d);
}

// user interface
void
dst3 (double *x, double *y, const struct aux *a) {
	s_dst3(x,y,1,a);
}

// strided one-dimensional DCT-4
inline static void
s_dct4_1d (double *x, double *y, int sy, const struct aux *a) {
	int n; // counter
	int N=a->N; // transform length
	double complex *t=a->t; // temporary buffer
	double complex *e=a->e; // exponent vector
	if (N==1) {
		// trivial case
		y[0] = 2*M_SQRT2*x[0];
		return;
	}
	// reduce to complex DFT of length N/2
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n)
		t[n] = 2**e++*(x[2*n]+I*x[N-1-2*n]);
	// do complex DFT in-place
	dft(t,t,a->sub1);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[sy*2*n] = 2*creal(*e++*t[n]);
		y[sy*(2*n+1)] = 2*creal(*e++*conj(t[N/2-1-n]));
	}
}

// strided DCT-4 of arbitrary dimension
inline static void
s_dct4 (double *x, double *y, int sy, const struct aux *a) {
	make_real_transform(x,y,sy,a,s_dct4_1d);
}

// user interface
void
dct4 (double *x, double *y, const struct aux *a) {
	s_dct4(x,y,1,a);
}

// strided one-dimensional DST-4
inline static void
s_dst4_1d (double *x, double *y, int sy, const struct aux *a) {
	int n; // counter
	int N=a->N; // transform length
	double complex *t=a->t; // temporary buffer
	double complex *e=a->e; // exponent vector
	if (N==1) {
		// trivial case
		y[0] = -2*M_SQRT2*x[0];
		return;
	}
	// reduce to complex DFT of length N/2
	// prepare sub-transform inputs
	for (n=0; n<N/2; ++n)
		t[n] = 2**e++*(x[2*n]-I*x[N-1-2*n]);
	// do complex DFT in-place
	dft(t,t,a->sub1);
	// recover results
	for (n=0; n<N/2; ++n) {
		y[sy*2*n] = 2*cimag(*e++*t[n]);
		y[sy*(2*n+1)] = 2*cimag(*e++*conj(t[N/2-1-n]));
	}
}

// strided DST-4 of arbitrary dimension
inline static void
s_dst4 (double *x, double *y, int sy, const struct aux *a) {
	make_real_transform(x,y,sy,a,s_dst4_1d);
}

// user interface
void
dst4 (double *x, double *y, const struct aux *a) {
	s_dst4(x,y,1,a);
}

// *** making of aux structures ***

// make aux structure for any transform of arbitrary dimension
// using its one-dimensional version
static struct aux *
mkaux_gen (int d, int *n, struct aux* (*mkaux_1d)(int N)) {
	struct aux *a;
	int p; // product of all transform lengths
	int i; // array index
	if (d==1)
		return (*mkaux_1d)(n[0]);
	else {
		p = 1;
		for (i=0; i<d; ++i)
			p *= n[i];
		a = malloc(sizeof(struct aux));
		a->N = p;
		a->t = malloc(p*sizeof(double complex));
		a->e = NULL;
		a->sub1 = mkaux_gen(d-1,n,mkaux_1d);
		a->sub2 = (*mkaux_1d)(n[d-1]);
		return a;
	}
}

// make aux for one-dimensional forward or inverse complex DFT
static struct aux *
mkaux_complex_1d (int N) {
	struct aux *a;
	int n;
	double complex *e;
	a = malloc(sizeof(struct aux));
	a->N = N;
	if (N>=16) {
		a->t = malloc(N*sizeof(double complex));
		a->e = malloc(N*sizeof(double complex));
		e = a->e;
		while (N>=16) {
			for (n=0; n<N/4; ++n) {
				*e++ = cexp(-2*M_PI*I*n/N);
				*e++ = cexp(-2*M_PI*I*3*n/N);
			}
			N /= 2;
		}
	} else {
		a->t = NULL;
		a->e = NULL;
	}
	a->sub1 = a->sub2 = NULL;
	return a;
}

// make aux for any-dimensional forward or inverse complex DFT
struct aux *
mkaux_complex (int d, int *n) {
	return mkaux_gen(d,n,mkaux_complex_1d);
}

// make aux for one-dimensional forward or inverse real DFT
struct aux *
mkaux_real_1d (int N) {
	struct aux *a;
	int n;
	double complex *e;
	a = malloc(sizeof(struct aux));
	a->N = N;
	a->t = NULL;
	if (N>=4) {
		a->e = malloc((N/4)*sizeof(double complex));
		e = a->e;
		for (n=0; n<N/4; ++n)
			*e++ = cexp(-2*M_PI*I*n/N);
		a->sub1 = mkaux_complex_1d(N/2);
	} else {
		a->e = NULL;
		a->sub1 = NULL;
	}
	a->sub2 = NULL;
	return a;
}

// make aux for one-dimensional Type-2 or Type-3 transforms
static struct aux *
mkaux_t2t3_1d (int N) {
	struct aux *a;
	int n;
	double complex *e;
	a = malloc(sizeof(struct aux));
	a->N = N;
	if (N>=2) {
		a->t = malloc(N*sizeof(double));
		a->e = malloc((N/2)*sizeof(double complex));
		e = a->e;
		for (n=0; n<N/2; ++n)
			*e++ = cexp(-2*M_PI*I*n/(4*N));
	} else {
		a->t = NULL;
		a->e = NULL;
	}
	a->sub1 = mkaux_real_1d(N);
	a->sub2 = NULL;
	return a;
}

// make aux for any-dimensional Type-2 or Type-3 transforms
struct aux *
mkaux_t2t3 (int d, int *n) {
	return mkaux_gen(d,n,mkaux_t2t3_1d);
}

// make aux for an one-dimensional Type-4 transform
static struct aux *
mkaux_t4_1d (int N) {
	struct aux *a;
	int n;
	double complex *e;
	a = malloc(sizeof(struct aux));
	a->N = N;
	if (N>=2) {
		a->t = malloc((N/2)*sizeof(double complex));
		a->e = malloc((N/2+N)*sizeof(double complex));
		e = a->e;
		for (n=0; n<N/2; ++n)
			*e++ = cexp(-2*M_PI*I*n/(2*N));
		for (n=0; n<N; ++n)
			*e++ = cexp(-2*M_PI*I*(2*n+1)/(8*N));
	} else {
		a->t = NULL;
		a->e = NULL;
	}
	a->sub1 = mkaux_complex_1d(N/2);
	a->sub2 = NULL;
	return a;
}

// make aux for an any-dimensional Type-4 transform
struct aux *
mkaux_t4 (int d, int *n) {
	return mkaux_gen(d,n,mkaux_t4_1d);
}

// free aux chain
void
free_aux (struct aux *a) {
	if (a->t)
		free(a->t);
	if (a->e)
		free(a->e);
	if (a->sub1)
		free_aux(a->sub1);
	if (a->sub2)
		free_aux(a->sub2);
	free(a);
}
