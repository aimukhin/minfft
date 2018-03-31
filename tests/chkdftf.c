#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#include "../dftf.h"

int
main (void) {
// test complex DFT for the given length
#if 0
	int n,N=8;
	complex float x[N],y[N];
	complex float *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// prepare exponent vector and do DFT
	e = mkexp_dftf(N);
	dftf(N,x,y,e);
	free(e);
	for (n=0; n<N; ++n)
		printf("%d %g %g\n",n,crealf(y[n]),cimagf(y[n]));
	// prepare exponent vector and do inverse DFT
	e = mkexp_idftf(N);
	idftf(N,y,x,e);
	free(e);
	for (n=0; n<N; ++n)
		printf("%d %g %g\n",n,crealf(x[n]),cimagf(x[n]));
#endif

// test real DFT for the given length
#if 0
	int n,N=8;
	float x[N],y[N];
	complex float *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// prepare exponent vector and do forward real DFT
	e = mkexp_realdftf(N);
	realdftf(N,x,y,e);
	free(e);
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,y[n]);
	// prepare exponent vector and do inverse real DFT
	e = mkexp_irealdftf(N);
	irealdftf(N,y,x,e);
	free(e);
	// print results
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,x[n]);
#endif

// test real symmetric DFTs for the given length
#if 0
	int n,N=4;
	float x[N],y[N];
	complex float *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// prepare exponent vector and do transform
//	e = mkexp_t2f(N); dct2f(N,x,y,e);
//	e = mkexp_t2f(N); dst2f(N,x,y,e);
//	e = mkexp_t3f(N); dct3f(N,x,y,e);
//	e = mkexp_t3f(N); dst3f(N,x,y,e);
//	e = mkexp_t4f(N); dct4f(N,x,y,e);
//	e = mkexp_t4f(N); dst4f(N,x,y,e);
	// print results
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,y[n]);
#endif

// inversion test for symmetric real DFTs of the given length
#if 0
	int n,N=4;
	float x[N],y[N];
	complex float *e;
	// init input vector
	for (n=0; n<N; ++n)
		x[n] = n+1;
	// do the transform and its inverse
//	e=mkexp_t2f(N); dct2f(N,x,y,e); free(e); e=mkexp_t3f(N); dct3f(N,y,x,e);
//	e=mkexp_t2f(N); dst2f(N,x,y,e); free(e); e=mkexp_t3f(N); dst3f(N,y,x,e);
//	e=mkexp_t4f(N); dct4f(N,x,y,e); dct4f(N,y,x,e);
//	e=mkexp_t4f(N); dst4f(N,x,y,e); dst4f(N,y,x,e);
	// print results
	for (n=0; n<N; ++n)
		printf("%d %g\n",n,x[n]);
#endif

// compare precision with plain complex DFT for the range of transform lengths
#if 0
#include <unistd.h>
	const int MAXBLK=65536*16;
	int n,N;
	float *x,*y;
	complex float *z,*w;
	complex float *e; // exponent vector
	double d,dmax; // maximum absolute error
	for (N=1; N<=MAXBLK; N*=2) {
//	N=1024; {
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (complex float*)malloc(8*N*sizeof(complex float));
		w = (complex float*)malloc(8*N*sizeof(complex float));
#if 0
		// real DFT
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_realdftf(N);
		realdftf(N,x,y,e);
		free(e);
		e = mkexp_dftf(N);
		dftf(N,z,w,e);
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
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=1; n<N/2; ++n) {
			z[n] = x[2*n]+I*x[2*n+1];
			z[N-n] = conjf(z[n]);
		}
		z[0] = x[0];
		if (N>1)
			z[N/2] = x[1];
		// do transforms
		e = mkexp_irealdftf(N);
		irealdftf(N,x,y,e);
		free(e);
		e = mkexp_idftf(N);
		idftf(N,z,w,e);
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
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<2*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = x[n];
		}
		// do transforms
		e = mkexp_t2f(N);
		dct2f(N,x,y,e);
		free(e);
		e = mkexp_dftf(4*N);
		dftf(4*N,z,w,e);
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
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<2*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = -x[n];
		}
		// do transforms
		e = mkexp_t2f(N);
		dst2f(N,x,y,e);
		free(e);
		e = mkexp_dftf(4*N);
		dftf(4*N,z,w,e);
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
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n) {
			z[n] = x[n];
			z[2*N-n] = -z[n];
		}
		z[N] = 0;
		for (n=0; n<2*N; ++n)
			z[2*N+n] = -z[n];
		// do transforms
		e = mkexp_t3f(N);
		dct3f(N,x,y,e);
		free(e);
		e = mkexp_dftf(4*N);
		dftf(4*N,z,w,e);
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
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n) {
			z[N-n] = x[n];
			z[N+n] = z[N-n];
		}
		z[0] = 0;
		for (n=0; n<2*N; ++n)
			z[2*N+n] = -z[n];
		// do transforms
		e = mkexp_t3f(N);
		dst3f(N,x,y,e);
		free(e);
		e = mkexp_dftf(4*N);
		dftf(4*N,z,w,e);
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
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<4*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = -z[2*n+1];
		}
		for (n=0; n<2*N; ++n)
			z[2*(4*N-1-n)+1] = z[2*n+1];
		// do transforms
		e = mkexp_t4f(N);
		dct4f(N,x,y,e);
		free(e);
		e = mkexp_dftf(8*N);
		dftf(8*N,z,w,e);
		free(e);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[2*n+1]));
			printf("%d %g\n",n,d);
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// DST-4
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<4*N; ++n)
			z[2*n] = 0;
		for (n=0; n<N; ++n) {
			z[2*n+1] = x[n];
			z[2*(2*N-1-n)+1] = z[2*n+1];
		}
		for (n=0; n<2*N; ++n)
			z[2*(4*N-1-n)+1] = -z[2*n+1];
		// do transforms
		e = mkexp_t4f(N);
		dst4f(N,x,y,e);
		free(e);
		e = mkexp_dftf(8*N);
		dftf(8*N,z,w,e);
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
	fftwf_plan p; // plan
	complex float *e; // exponent vector
	for (N=1; N<=MAXBLK; N*=2) {
//	N=1024; {
#if 0
		// complex DFT
		complex float *x,*y,*z,*w;
		x = (complex float*)malloc(N*sizeof(complex float));
		y = (complex float*)malloc(N*sizeof(complex float));
		z = (complex float*)malloc(N*sizeof(complex float));
		w = (complex float*)malloc(N*sizeof(complex float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = \
			(float)rand()/RAND_MAX-0.5+I*((float)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_dftf(N);
		dftf(N,x,y,e);
		p = fftwf_plan_dft_1d(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT); 
		fftwf_execute(p);
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
		complex float *x,*y,*z,*w;
		x = (complex float*)malloc(N*sizeof(complex float));
		y = (complex float*)malloc(N*sizeof(complex float));
		z = (complex float*)malloc(N*sizeof(complex float));
		w = (complex float*)malloc(N*sizeof(complex float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = \
			(float)rand()/RAND_MAX-0.5+I*((float)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_idftf(N);
		idftf(N,x,y,e);
		p = fftwf_plan_dft_1d(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT); 
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_realdftf(N);
		realdftf(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_R2HC,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		z[0] = x[0];
		if (N>1)
			z[N/2] = x[1];
		for (n=1; n<N/2; ++n) {
			z[n] = x[2*n];
			z[N-n] = x[2*n+1];
		}
		// do transforms
		e = mkexp_irealdftf(N);
		irealdftf(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_HC2R,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t2f(N);
		dct2f(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_REDFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t2f(N);
		dst2f(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_RODFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t3f(N);
		dct3f(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_REDFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[N-1-n];
		// do transforms
		e = mkexp_t3f(N);
		dst3f(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_RODFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t4f(N);
		dct4f(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		float *x,*y,*z,*w;
		x = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		z = (float*)malloc(N*sizeof(float));
		w = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		for (n=0; n<N; ++n)
			z[n] = x[n];
		// do transforms
		e = mkexp_t4f(N);
		dst4f(N,x,y,e);
		p = fftwf_plan_r2r_1d(N,z,w,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
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
		fftwf_destroy_plan(p);
		printf("%d\tmax = %g\n",N,dmax);
	}
#endif

// performance test of complex transforms
#if 0
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	complex float *z = (complex float*)malloc(MAXBLK*sizeof(complex float));
	complex float *w = (complex float*)malloc(MAXBLK*sizeof(complex float));
	complex float *e;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare test vector
		for (n=0; n<N; ++n)
			z[n] = \
			(float)rand()/RAND_MAX-0.5+I*((float)rand()/RAND_MAX-0.5);
		// prepare exponent vector
		e = mkexp_dftf(N);
//		e = mkexp_idftf(N);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				dftf(N,z,w,e);
//				idftf(N,z,w,e);
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
	fftwf_complex *z,*w;
	fftwf_plan p;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	z = (fftwf_complex*)fftw_malloc(MAXBLK*sizeof(fftwf_complex));
	w = (fftwf_complex*)fftw_malloc(MAXBLK*sizeof(fftwf_complex));
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare plan
		p = fftwf_plan_dft_1d(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
//		p = fftwf_plan_dft_1d(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// prepare test vector
		for (n=0; n<N; ++n)
			z[n] = \
			(float)rand()/RAND_MAX-0.5+I*((float)rand()/RAND_MAX-0.5);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				fftwf_execute(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%10d %10d %10.4g %10.4g\t\t%g\n",N,R*T,avg,stdd,cabs(w[0]));
		fftwf_destroy_plan(p);
	}
#endif

// performance test of real transforms
#if 0
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	float *x = (float*)malloc(MAXBLK*sizeof(float));
	float *y = (float*)malloc(MAXBLK*sizeof(float));
	complex float *e;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX-0.5;
		// prepare exponent vector
//		e = mkexp_realdftf(N);
//		e = mkexp_irealdftf(N);
//		e = mkexp_t2f(N);
//		e = mkexp_t3f(N);
//		e = mkexp_t4f(N);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				realdftf(N,x,y,e);
//				irealdftf(N,x,y,e);
//				dct2f(N,x,y,e);
//				dst2f(N,x,y,e);
//				dct3f(N,x,y,e);
//				dst3f(N,x,y,e);
//				dct4f(N,x,y,e);
//				dst4f(N,x,y,e);
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
	float *x,*y;
	fftwf_plan p;
	int N,n,r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	x = (float*)fftw_malloc(MAXBLK*sizeof(float));
	y = (float*)fftw_malloc(MAXBLK*sizeof(float));
	for (N=1; N<=MAXBLK; N*=2) {
		// prepare plan
		// realdft
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_R2HC,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// inverse realdft
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_HC2R,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dct-2
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dst-2
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dct-3
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_REDFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dst-3
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_RODFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dct-4
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// dst-4
//		p = fftwf_plan_r2r_1d(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (float)rand()/RAND_MAX;
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				fftwf_execute(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%10d %10d %10.4g %10.4g\t\t%g\n",N,R*T,avg,stdd,y[0]);
		fftwf_destroy_plan(p);
	}
#endif

// compare precision of forward and inverse transforms
#if 0
#include <unistd.h>
	const int MAXBLK=65536;
	int n,N;
	double d,dmax; // maximum absolute error
	complex float *e; // exponent vector
	for (N=1; N<=MAXBLK; N*=2) {
#if 0
		// complex forward and inverse DFT
		complex float *x,*x0,*y;
		x = (complex float*)malloc(N*sizeof(complex float));
		x0 = (complex float*)malloc(N*sizeof(complex float));
		y = (complex float*)malloc(N*sizeof(complex float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = \
			(float)rand()/RAND_MAX-0.5+I*((float)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_dftf(N);
		dftf(N,x,y,e);
		free(e);
		e = mkexp_idftf(N);
		idftf(N,y,x,e);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_realdftf(N);
		realdftf(N,x,y,e);
		free(e);
		e = mkexp_irealdftf(N);
		irealdftf(N,y,x,e);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t2f(N);
		dct2f(N,x,y,e);
		free(e);
		e = mkexp_t3f(N);
		dct3f(N,y,x,e);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t2f(N);
		dst2f(N,x,y,e);
		free(e);
		e = mkexp_t3f(N);
		dst3f(N,y,x,e);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t4f(N);
		dct4f(N,x,y,e);
		dct4f(N,y,x,e);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		e = mkexp_t4f(N);
		dst4f(N,x,y,e);
		dst4f(N,y,x,e);
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
	fftwf_plan p;
	for (N=1; N<=MAXBLK; N*=2) {
#if 0
		// complex forward and inverse DFT
		complex float*x,*x0,*y;
		x = (complex float*)malloc(N*sizeof(complex float));
		x0 = (complex float*)malloc(N*sizeof(complex float));
		y = (complex float*)malloc(N*sizeof(complex float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = \
			(float)rand()/RAND_MAX-0.5+I*((float)rand()/RAND_MAX-0.5);
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftwf_plan_dft_1d(N,x,y,FFTW_FORWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_dft_1d(N,y,x,FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftwf_plan_r2r_1d(N,x,y,FFTW_R2HC,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(N,y,x,FFTW_HC2R,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftwf_plan_r2r_1d(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(N,y,x,FFTW_REDFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftwf_plan_r2r_1d(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(N,y,x,FFTW_RODFT01,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftwf_plan_r2r_1d(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(N,y,x,FFTW_REDFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
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
		float *x,*x0,*y;
		x = (float*)malloc(N*sizeof(float));
		x0 = (float*)malloc(N*sizeof(float));
		y = (float*)malloc(N*sizeof(float));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x0[n] = (float)rand()/RAND_MAX;
		for (n=0; n<N; ++n)
			x[n] = x0[n];
		// do transforms
		p = fftwf_plan_r2r_1d(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
		p = fftwf_plan_r2r_1d(N,y,x,FFTW_RODFT11,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
		fftwf_execute(p);
		fftwf_destroy_plan(p);
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
