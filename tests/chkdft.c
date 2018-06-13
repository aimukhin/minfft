#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#include "../dft.h"

int
main (void) {
// test complex DFT for the given length
#if 0
	int n,N=8;
	double complex x[N],y[N];
	double complex *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = (n+1)*(1+I);
	// prepare exponent vector and do DFT
	e = mkexp_dft(N);
	dft(N,x,y,e);
	free(e);
	for (n=0; n<N; ++n)
		printf("%d %g %g\n",n,creal(y[n]),cimag(y[n]));
	// prepare exponent vector and do inverse DFT
	e = mkexp_idft(N);
	idft(N,y,x,e);
	free(e);
	for (n=0; n<N; ++n)
		printf("%d %g %g\n",n,creal(x[n]),cimag(x[n]));
#endif

// test real DFT for the given length
#if 0
	int n,N=8;
	double x[N],y[N];
	double complex *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// prepare exponent vector and do forward real DFT
	e = mkexp_realdft(N);
	realdft(N,x,y,e);
	free(e);
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,y[n]);
	// prepare exponent vector and do inverse real DFT
	e = mkexp_irealdft(N);
	irealdft(N,y,x,e);
	free(e);
	// print results
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,x[n]);
#endif

// test real symmetric DFTs for the given length
#if 0
	int n,N=4;
	double x[N],y[N];
	double complex *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// prepare exponent vector and do transform
//	e = mkexp_t2(N); dct2(N,x,y,e);
//	e = mkexp_t2(N); dst2(N,x,y,e);
//	e = mkexp_t3(N); dct3(N,x,y,e);
//	e = mkexp_t3(N); dst3(N,x,y,e);
//	e = mkexp_t4(N); dct4(N,x,y,e);
//	e = mkexp_t4(N); dst4(N,x,y,e);
	// print results
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,y[n]);
#endif

// inversion test for symmetric real DFTs of the given length
#if 0
	int n,N=4;
	double x[N],y[N];
	double complex *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// do the transform and its inverse
//	e=mkexp_t2(N); dct2(N,x,y,e); free(e); e=mkexp_t3(N); dct3(N,y,x,e);
//	e=mkexp_t2(N); dst2(N,x,y,e); free(e); e=mkexp_t3(N); dst3(N,y,x,e);
//	e=mkexp_t4(N); dct4(N,x,y,e); dct4(N,y,x,e);
//	e=mkexp_t4(N); dst4(N,x,y,e); dst4(N,y,x,e);
	// print results
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,x[n]);
#endif

// compare precision with plain complex DFT for the range of transform lengths
#if 0
#include <unistd.h>
	const int MAXBLK=65536*16;
	int n,N;
	double *x,*y;
	double complex *z,*w;
	double complex *e; // exponent vector
	double d,dmax; // maximum absolute error
	for (N=1; N<=MAXBLK; N*=2) {
//	N=1024; {
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double complex*)malloc(8*N*sizeof(double complex));
		w = (double complex*)malloc(8*N*sizeof(double complex));
#if 0
		// real DFT
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_realdft(N);
		realdft(N,x,y,e);
		free(e);
		e = mkexp_dft(N);
		dft(N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=1; n<N/2; ++n) {
			d = log10(cabs(y[2*n]+I*y[2*n+1]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
		// boundary cases
		d = log10(cabs(y[0]-w[0]));
		dmax = (d>dmax)?d:dmax;
		if (N>1) {
			d = log10(cabs(y[1]-w[N/2]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// inverse real DFT
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=1; n<N/2; ++n) {
			z[n] = x[2*n]+I*x[2*n+1];
			z[N-n] = conj(z[n]);
		}
		z[0] = x[0];
		if (N>1)
			z[N/2] = x[1];
		// do transforms
		e = mkexp_irealdft(N);
		irealdft(N,x,y,e);
		free(e);
		e = mkexp_idft(N);
		idft(N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-2
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<2*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = x[n];
		}
		// do transforms
		e = mkexp_t2(N);
		dct2(N,x,y,e);
		free(e);
		e = mkexp_dft(4*N);
		dft(4*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-2
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<2*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = -x[n];
		}
		// do transforms
		e = mkexp_t2(N);
		dst2(N,x,y,e);
		free(e);
		e = mkexp_dft(4*N);
		dft(4*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(I*y[n]-w[N-n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-3
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n) {
			z[n] = x[n];
			z[2*N-n] = -z[n];
		}
		z[N] = 0;
		for (n=0; n<2*N; ++n)
			z[2*N+n] = -z[n];
		// do transforms
		e = mkexp_t3(N);
		dct3(N,x,y,e);
		free(e);
		e = mkexp_dft(4*N);
		dft(4*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[2*n+1]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-3
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n) {
			z[N-n] = x[n];
			z[N+n] = z[N-n];
		}
		z[0] = 0;
		for (n=0; n<2*N; ++n)
			z[2*N+n] = -z[n];
		// do transforms
		e = mkexp_t3(N);
		dst3(N,x,y,e);
		free(e);
		e = mkexp_dft(4*N);
		dft(4*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(I*y[n]-w[2*n+1]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-4
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<4*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = -z[2*n+1];
		}
		for (n=0; n<2*N; ++n)
			z[2*(4*N-1-n)+1] = z[2*n+1];
		// do transforms
		e = mkexp_t4(N);
		dct4(N,x,y,e);
		free(e);
		e = mkexp_dft(8*N);
		dft(8*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[2*n+1]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-4
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<4*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = z[2*n+1];
		}
		for (n=0; n<2*N; ++n)
			z[2*(4*N-1-n)+1] = -z[2*n+1];
		// do transforms
		e = mkexp_t4(N);
		dst4(N,x,y,e);
		free(e);
		e = mkexp_dft(8*N);
		dft(8*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(I*y[n]-w[2*n+1]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(y);
		free(z);
		free(w);
		printf("%d\tmax = %g\n",N,dmax);
	}
#endif

// compare precision with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXBLK=65536*16;
	int n,N;
	double d,dmax; // maximum absolute error
	fftw_plan p; // plan
	double complex *e; // exponent vector
	for (N=1; N<=MAXBLK; N*=2) {
//	N=1024; {
#if 0
		// complex DFT
		double complex *x,*y,*z,*w;
		x = (double complex*)malloc(N*sizeof(double complex));
		y = (double complex*)malloc(N*sizeof(double complex));
		z = (double complex*)malloc(N*sizeof(double complex));
		w = (double complex*)malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+I*((double)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_dft(N);
		dft(N,x,y,e);
		p = fftw_plan_dft_1d(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y,*z,*w;
		x = (double complex*)malloc(N*sizeof(double complex));
		y = (double complex*)malloc(N*sizeof(double complex));
		z = (double complex*)malloc(N*sizeof(double complex));
		w = (double complex*)malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+I*((double)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_idft(N);
		idft(N,x,y,e);
		p = fftw_plan_dft_1d(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// real DFT
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_realdft(N);
		realdft(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_R2HC,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N/2; ++n) {
			d = log10(fabs(y[2*n]-w[n]));
//			printf("%d %g\n",2*n,d);
			dmax = (d>dmax)?d:dmax;
		}
		if (N>1) {
			d = log10(fabs(y[1]-w[N/2]));
//			printf("%d %g\n",1,d);
			dmax = (d>dmax)?d:dmax;
		}
		for (n=1; n<N/2; ++n) {
			d = log10(fabs(y[2*n+1]-w[N-n]));
//			printf("%d %g\n",2*n+1,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// inverse real DFT
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		z[0] = x[0];
		if (N>1)
			z[N/2] = x[1];
		for (n=1; n<N/2; ++n) {
			z[n] = x[2*n];
			z[N-n] = x[2*n+1];
		}
		// do transforms
		e = mkexp_irealdft(N);
		irealdft(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_HC2R,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-2
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t2(N);
		dct2(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-2
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t2(N);
		dst2(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(-y[N-1-n]-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-3
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t3(N);
		dct3(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]/2-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-3
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[N-1-n];
		// do transforms
		e = mkexp_t3(N);
		dst3(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(-y[n]/2-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-4
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t4(N);
		dct4(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]/2-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-4
		double *x,*y,*z,*w;
		x = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		z = (double*)malloc(N*sizeof(double));
		w = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t4(N);
		dst4(N,x,y,e);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(-y[n]/2-w[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(y);
		free(z);
		free(w);
		free(e);
		fftw_destroy_plan(p);
		printf("%d\tmax = %g\n",N,dmax);
	}
#endif

// performance test of complex transforms
#if 0
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	double complex *z = (double complex*)malloc(MAXBLK*sizeof(double complex));
	double complex *w = (double complex*)malloc(MAXBLK*sizeof(double complex));
	double complex *e;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare test vector
		for (n=0; n<N; ++n)
			z[n] = \
			(double)rand()/RAND_MAX-0.5+I*((double)rand()/RAND_MAX-0.5);
		// prepare exponent vector
		e = mkexp_dft(N);
//		e = mkexp_idft(N);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				dft(N,z,w,e);
//				idft(N,z,w,e);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%10d %10d %10.4g %10.4g\t\t%g\n",N,R*T,avg,stdd,cabs(w[0]));
		free(e);
	}
#endif

// performance test of complex FFTW
#if 0
#include <fftw3.h>
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	fftw_complex *z,*w;
	fftw_plan p;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	z = (fftw_complex*)fftw_malloc(MAXBLK*sizeof(fftw_complex));
	w = (fftw_complex*)fftw_malloc(MAXBLK*sizeof(fftw_complex));
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare plan
		p = fftw_plan_dft_1d(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
//		p = fftw_plan_dft_1d(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// prepare test vector
		for (n=0; n<N; ++n)
			z[n] = \
			(double)rand()/RAND_MAX-0.5+I*((double)rand()/RAND_MAX-0.5);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				fftw_execute(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%10d %10d %10.4g %10.4g\t\t%g\n",N,R*T,avg,stdd,cabs(w[0]));
		fftw_destroy_plan(p);
	}
#endif

// performance test of real transforms
#if 0
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	double *x = (double*)malloc(MAXBLK*sizeof(double));
	double *y = (double*)malloc(MAXBLK*sizeof(double));
	double complex *e;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare exponent vector
//		e = mkexp_realdft(N);
//		e = mkexp_irealdft(N);
//		e = mkexp_t2(N);
//		e = mkexp_t3(N);
//		e = mkexp_t4(N);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				realdft(N,x,y,e);
//				irealdft(N,x,y,e);
//				dct2(N,x,y,e);
//				dst2(N,x,y,e);
//				dct3(N,x,y,e);
//				dst3(N,x,y,e);
//				dct4(N,x,y,e);
//				dst4(N,x,y,e);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%10d %10d %10.4g %10.4g\t\t%g\n",N,R*T,avg,stdd,y[0]);
		free(e);
	}
#endif

// performance test of real FFTW
#if 0
#include <fftw3.h>
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	double *x,*y;
	fftw_plan p;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	x = (double*)fftw_malloc(MAXBLK*sizeof(double));
	y = (double*)fftw_malloc(MAXBLK*sizeof(double));
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare plan
		// realdft
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_R2HC,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// inverse realdft
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_HC2R,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dct-2
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dst-2
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dct-3
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dst-3
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dct-4
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dst-4
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX;
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				fftw_execute(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%10d %10d %10.4g %10.4g\t\t%g\n",N,R*T,avg,stdd,y[0]);
		fftw_destroy_plan(p);
	}
#endif

// compare precision of forward and inverse transforms
#if 0
#include <unistd.h>
	const int MAXBLK=65536;
	int n,N;
	double d,dmax; // maximum absolute error
	double complex *e; // exponent vector
	for (N=1; N<=MAXBLK; N*=2) {
#if 0
		// complex forward and inverse DFT
		double complex *x,*x0,*y;
		x = (double complex*)malloc(N*sizeof(double complex));
		x0 = (double complex*)malloc(N*sizeof(double complex));
		y = (double complex*)malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = \
			(double)rand()/RAND_MAX-0.5+I*((double)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_dft(N);
		dft(N,x,y,e);
		free(e);
		e = mkexp_idft(N);
		idft(N,y,x,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(x[n]-N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// real forward and inverse DFT
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_realdft(N);
		realdft(N,x,y,e);
		free(e);
		e = mkexp_irealdft(N);
		irealdft(N,y,x,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t2(N);
		dct2(N,x,y,e);
		free(e);
		e = mkexp_t3(N);
		dct3(N,y,x,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-4*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t2(N);
		dst2(N,x,y,e);
		free(e);
		e = mkexp_t3(N);
		dst3(N,y,x,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-4*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-4
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t4(N);
		dct4(N,x,y,e);
		dct4(N,y,x,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-8*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-4
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t4(N);
		dst4(N,x,y,e);
		dst4(N,y,x,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-8*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(x0);
		free(y);
		printf("%d\tmax = %g\n",N,dmax);
	}
#endif

// compare precision of forward and inverse transforms for FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXBLK=65536;
	int n,N;
	double d,dmax; // maximum absolute error
	fftw_plan p;
	for (N=1; N<=MAXBLK; N*=2) {
#if 0
		// complex forward and inverse DFT
		double complex *x,*x0,*y;
		x = (double complex*)malloc(N*sizeof(double complex));
		x0 = (double complex*)malloc(N*sizeof(double complex));
		y = (double complex*)malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = \
			(double)rand()/RAND_MAX-0.5+I*((double)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftw_plan_dft_1d(N,x,y,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_dft_1d(N,y,x,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(x[n]-N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// real forward and inverse DFT
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftw_plan_r2r_1d(N,x,y,FFTW_R2HC,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_1d(N,y,x,FFTW_HC2R,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_1d(N,y,x,FFTW_REDFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-2*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_1d(N,y,x,FFTW_RODFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-2*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DCT-4
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_1d(N,y,x,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-2*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-4
		double *x,*x0,*y;
		x = (double*)malloc(N*sizeof(double));
		x0 = (double*)malloc(N*sizeof(double));
		y = (double*)malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (double)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		p = fftw_plan_r2r_1d(N,y,x,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftw_execute(p);
		fftw_destroy_plan(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(x[n]-2*N*x0[n]));
//			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(x0);
		free(y);
		printf("%d\tmax = %g\n",N,dmax);
	}
#endif
}
