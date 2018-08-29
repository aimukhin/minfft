#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#include "../minfft.h"

int
main (void) {

// one-dimensional DFTs

// compare results of one-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXN=65536*16;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	fftw_plan p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex DFT
		double complex *x,*y,*z,*w;
		x = malloc(N*sizeof(double complex));
		y = malloc(N*sizeof(double complex));
		z = malloc(N*sizeof(double complex));
		w = malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		p = fftw_plan_dft_1d(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y,*z,*w;
		x = malloc(N*sizeof(double complex));
		y = malloc(N*sizeof(double complex));
		z = malloc(N*sizeof(double complex));
		w = malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_invdft(x,y,a);
		p = fftw_plan_dft_1d(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
#endif
#if 0
		// DCT-2
		double *x,*y,*z,*w;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		w = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dct2(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT10,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-2
		double *x,*y,*z,*w;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		w = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dst2(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT10,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DCT-3
		double *x,*y,*z,*w;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		w = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dct3(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT01,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-3
		double *x,*y,*z,*w;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		w = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dst3(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT01,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DCT-4
		double *x,*y,*z,*w;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		w = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dct4(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT11,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-4
		double *x,*y,*z,*w;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		w = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dst4(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT11,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = cabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%8d %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		fftw_destroy_plan(p);
	}
#endif

// compare results of one-dimensional transforms with Kiss FFT
#if 0
#include <unistd.h>
#include "kiss_fft.h"
	const int MAXN=65536*16;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	kiss_fft_cfg cfg; // Kiss FFT config
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex DFT
		double complex *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*sizeof(double complex));
		y = malloc(N*sizeof(double complex));
		z = malloc(N*sizeof(kiss_fft_cpx));
		w = malloc(N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n) {
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		cfg = kiss_fft_alloc(N,0,NULL,NULL);
		kiss_fft(cfg,z,w);
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*sizeof(double complex));
		y = malloc(N*sizeof(double complex));
		z = malloc(N*sizeof(kiss_fft_cpx));
		w = malloc(N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n) {
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_invdft(x,y,a);
		cfg = kiss_fft_alloc(N,1,NULL,NULL);
		kiss_fft(cfg,z,w);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = cabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%8d %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		free(cfg);
	}
#endif

// compare forward and inverse one-dimensional transforms
// with the identity operator
#if 0
#include <unistd.h>
	const int MAXN=65536*16;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// forward and inverse complex DFT
		double complex *x,*y,*z;
		x = malloc(N*sizeof(double complex));
		y = malloc(N*sizeof(double complex));
		z = malloc(N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		minfft_invdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(x[n]-z[n]/N);
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		double *x,*y,*z;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dct2(x,y,a);
		minfft_dct3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		double *x,*y,*z;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dst2(x,y,a);
		minfft_dst3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-4
		double *x,*y,*z;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dct4(x,y,a);
		minfft_dct4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-4
		double *x,*y,*z;
		x = malloc(N*sizeof(double));
		y = malloc(N*sizeof(double));
		z = malloc(N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dst4(x,y,a);
		minfft_dst4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%8d %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		minfft_free_aux(a);
	}
#endif

// performance test of one-dimensional transforms
#if 0
	const int MINT=1;
	const int MAXN=65536*16;
	const int R=10;
	minfft_aux *a;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		double complex *x = malloc(N*sizeof(double complex));
		double complex *y = malloc(N*sizeof(double complex));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// prepare aux data
		a = minfft_mkaux_dft_1d(N);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				minfft_dft(x,y,a);
//				minfft_invdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = malloc(N*sizeof(double));
		double *y = malloc(N*sizeof(double));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare aux data
//		a = minfft_mkaux_t2t3_1d(N);
//		a = minfft_mkaux_t4_1d(N);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				minfft_dct2(x,y,a);
//				minfft_dst2(x,y,a);
//				minfft_dct3(x,y,a);
//				minfft_dst3(x,y,a);
//				minfft_dct4(x,y,a);
//				minfft_dst4(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		minfft_free_aux(a);
	}
#endif

// performance test of one-dimensional FFTW
#if 0
#include <fftw3.h>
	const int MINT=1;
	const int MAXN=65536*16;
	const int R=10;
	fftw_plan p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		fftw_complex *x = fftw_malloc(N*sizeof(double complex));
		fftw_complex *y = fftw_malloc(N*sizeof(double complex));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// prepare plan
//		p = fftw_plan_dft_1d(N,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
//		p = fftw_plan_dft_1d(N,x,y,FFTW_BACKWARD,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
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
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = fftw_malloc(N*sizeof(double));
		double *y = fftw_malloc(N*sizeof(double));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
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
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		fftw_destroy_plan(p);
	}
#endif

// performance test of one-dimensional Kiss FFT
#if 0
#include "kiss_fft.h"
	const int MINT=1;
	const int MAXN=65536*16;
	const int R=10;
	kiss_fft_cfg cfg;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
		// complex transforms
		kiss_fft_cpx *x = malloc(N*sizeof(kiss_fft_cpx));
		kiss_fft_cpx *y = malloc(N*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N; ++n) {
			x[n].r = (double)rand()/RAND_MAX-0.5;
			x[n].i = (double)rand()/RAND_MAX-0.5;
		}
		// prepare config
//		cfg = kiss_fft_alloc(N,0,NULL,NULL);
//		cfg = kiss_fft_alloc(N,1,NULL,NULL);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fft(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,hypot(y[0].r,y[0].i));
		free(x);
		free(y);
		free(cfg);
#endif

// two-dimensional DFTs

// compare results of two-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXN=1024;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	fftw_plan p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex DFT
		double complex *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double complex));
		y = malloc(N*N*sizeof(double complex));
		z = malloc(N*N*sizeof(double complex));
		w = malloc(N*N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_dft(x,y,a);
		p = fftw_plan_dft_2d(N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double complex));
		y = malloc(N*N*sizeof(double complex));
		z = malloc(N*N*sizeof(double complex));
		w = malloc(N*N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_invdft(x,y,a);
		p = fftw_plan_dft_2d(N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
#endif
#if 0
		// DCT-2
		double *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		w = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dct2(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-2
		double *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		w = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dst2(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DCT-3
		double *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		w = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dct3(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-3
		double *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		w = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dst3(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DCT-4
		double *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		w = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dct4(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-4
		double *x,*y,*z,*w;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		w = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dst4(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = cabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = cabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%5d**2 %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		fftw_destroy_plan(p);
	}
#endif

// compare results of two-dimensional transforms with Kiss FFT
#if 0
#include <unistd.h>
#include "kiss_fftnd.h"
	const int MAXN=1024;
	int N,n,Ns[2];
	double d,dmax,v,vmax; // current and maximum absolute values
	kiss_fftnd_cfg cfg; // Kiss FFT config
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex DFT
		double complex *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*sizeof(double complex));
		y = malloc(N*N*sizeof(double complex));
		z = malloc(N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n) {
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_dft(x,y,a);
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftnd_alloc(Ns,2,0,NULL,NULL);
		kiss_fftnd(cfg,z,w);
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*sizeof(double complex));
		y = malloc(N*N*sizeof(double complex));
		z = malloc(N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n) {
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_invdft(x,y,a);
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftnd_alloc(Ns,2,1,NULL,NULL);
		kiss_fftnd(cfg,z,w);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = cabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = cabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%5d**2 %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		free(cfg);
	}
#endif

// compare forward and inverse two-dimensional transforms
// with the identity operator
#if 0
#include <unistd.h>
	const int MAXN=1024;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// forward and inverse complex DFT
		double complex *x,*y,*z;
		x = malloc(N*N*sizeof(double complex));
		y = malloc(N*N*sizeof(double complex));
		z = malloc(N*N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_dft(x,y,a);
		minfft_invdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = cabs(x[n]-z[n]/(N*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		double *x,*y,*z;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dct2(x,y,a);
		minfft_dct3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		double *x,*y,*z;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dst2(x,y,a);
		minfft_dst3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-4
		double *x,*y,*z;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dct4(x,y,a);
		minfft_dct4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-4
		double *x,*y,*z;
		x = malloc(N*N*sizeof(double));
		y = malloc(N*N*sizeof(double));
		z = malloc(N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dst4(x,y,a);
		minfft_dst4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%5d**2 %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		minfft_free_aux(a);
	}
#endif

// performance test of two-dimensional transforms
#if 0
	const int MINT=1;
	const int MAXN=1024;
	const int R=10;
	minfft_aux *a;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		double complex *x = malloc(N*N*sizeof(double complex));
		double complex *y = malloc(N*N*sizeof(double complex));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// prepare aux data
		a = minfft_mkaux_dft_2d(N,N);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				minfft_dft(x,y,a);
//				minfft_invdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = malloc(N*N*sizeof(double));
		double *y = malloc(N*N*sizeof(double));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] =(double)rand()/RAND_MAX-0.5;
		// prepare aux data
//		a = minfft_mkaux_t2t3_2d(N,N);
//		a = minfft_mkaux_t4_2d(N,N);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				minfft_dct2(x,y,a);
//				minfft_dst2(x,y,a);
//				minfft_dct3(x,y,a);
//				minfft_dst3(x,y,a);
//				minfft_dct4(x,y,a);
//				minfft_dst4(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		minfft_free_aux(a);
	}
#endif

// performance test of two-dimensional FFTW
#if 0
#include <fftw3.h>
	const int MINT=1;
	const int MAXN=1024;
	const int R=10;
	fftw_plan p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		fftw_complex *x = fftw_malloc(N*N*sizeof(fftw_complex));
		fftw_complex *y = fftw_malloc(N*N*sizeof(fftw_complex));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// prepare plan
//		p = fftw_plan_dft_2d(N,N,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
//		p = fftw_plan_dft_2d(N,N,x,y,FFTW_BACKWARD,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
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
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = fftw_malloc(N*N*sizeof(double));
		double *y = fftw_malloc(N*N*sizeof(double));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = fftw_plan_r2r_2d(N,N,x,y,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_2d(N,N,x,y,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_2d(N,N,x,y,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_2d(N,N,x,y,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_2d(N,N,x,y,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_2d(N,N,x,y,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
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
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		fftw_destroy_plan(p);
	}
#endif

// performance test of two-dimensional Kiss FFT
#if 0
#include "kiss_fftnd.h"
	const int MINT=1;
	const int MAXN=1024;
	const int R=10;
	kiss_fftnd_cfg cfg;
	int N,n,Ns[2];
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
		// complex transforms
		kiss_fft_cpx *x = malloc(N*N*sizeof(kiss_fft_cpx));
		kiss_fft_cpx *y = malloc(N*N*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N*N; ++n) {
			x[n].r = (double)rand()/RAND_MAX-0.5;
			x[n].i = (double)rand()/RAND_MAX-0.5;
		}
		// prepare config
		Ns[0] = Ns[1] = N;
//		cfg = kiss_fftnd_alloc(Ns,2,0,NULL,NULL);
//		cfg = kiss_fftnd_alloc(Ns,2,1,NULL,NULL);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftnd(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,hypot(y[0].r,y[0].i));
		free(x);
		free(y);
		free(cfg);
	}
#endif

// three-dimensional DFTs

// compare results of three-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXN=128;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	fftw_plan p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex DFT
		double complex *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double complex));
		y = malloc(N*N*N*sizeof(double complex));
		z = malloc(N*N*N*sizeof(double complex));
		w = malloc(N*N*N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_dft(x,y,a);
		p = fftw_plan_dft_3d(N,N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double complex));
		y = malloc(N*N*N*sizeof(double complex));
		z = malloc(N*N*N*sizeof(double complex));
		w = malloc(N*N*N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_invdft(x,y,a);
		p = fftw_plan_dft_3d(N,N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
#endif
#if 0
		// DCT-2
		double *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		w = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dct2(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-2
		double *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		w = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dst2(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DCT-3
		double *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		w = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dct3(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-3
		double *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		w = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dst3(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DCT-4
		double *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		w = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dct4(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
#if 0
		// DST-4
		double *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		w = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dst4(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		fftw_execute(p);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = cabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%5d**3 %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
	}
#endif

// compare results of three-dimensional transforms with Kiss FFT
#if 0
#include <unistd.h>
#include "kiss_fftnd.h"
	const int MAXN=128;
	int N,n,Ns[3];
	double d,dmax,v,vmax; // current and maximum absolute values
	kiss_fftnd_cfg cfg; // Kiss FFT config
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex DFT
		double complex *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*N*sizeof(double complex));
		y = malloc(N*N*N*sizeof(double complex));
		z = malloc(N*N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n) {
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_dft(x,y,a);
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftnd_alloc(Ns,3,0,NULL,NULL);
		kiss_fftnd(cfg,z,w);
#endif
#if 0
		// inverse complex DFT
		double complex *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*N*sizeof(double complex));
		y = malloc(N*N*N*sizeof(double complex));
		z = malloc(N*N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n) {
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_invdft(x,y,a);
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftnd_alloc(Ns,3,1,NULL,NULL);
		kiss_fftnd(cfg,z,w);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = cabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%5d**3 %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		free(cfg);
	}
#endif

// compare forward and inverse three-dimensional transforms
// with the identity operator
#if 0
#include <unistd.h>
	const int MAXN=128;
	int N,n;
	double d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// forward and inverse complex DFT
		double complex *x,*y,*z;
		x = malloc(N*N*N*sizeof(double complex));
		y = malloc(N*N*N*sizeof(double complex));
		z = malloc(N*N*N*sizeof(double complex));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_dft(x,y,a);
		minfft_invdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(x[n]-z[n]/(N*N*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		double *x,*y,*z;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dct2(x,y,a);
		minfft_dct3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		double *x,*y,*z;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dst2(x,y,a);
		minfft_dst3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-4
		double *x,*y,*z;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dct4(x,y,a);
		minfft_dct4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-4
		double *x,*y,*z;
		x = malloc(N*N*N*sizeof(double));
		y = malloc(N*N*N*sizeof(double));
		z = malloc(N*N*N*sizeof(double));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dst4(x,y,a);
		minfft_dst4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = cabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = cabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%5d**3 %g\n",N,dmax/vmax);
		// free resources
		free(x);
		free(y);
		free(z);
		minfft_free_aux(a);
	}
#endif

// performance test of three-dimensional transforms
#if 0
	const int MINT=1;
	const int MAXN=128;
	const int R=10;
	minfft_aux *a;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		double complex *x = malloc(N*N*N*sizeof(double complex));
		double complex *y = malloc(N*N*N*sizeof(double complex));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// prepare aux data
		a = minfft_mkaux_dft_3d(N,N,N);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				minfft_dft(x,y,a);
//				minfft_invdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = malloc(N*N*N*sizeof(double));
		double *y = malloc(N*N*N*sizeof(double));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] =(double)rand()/RAND_MAX-0.5;
		// prepare aux data
//		a = minfft_mkaux_t2t3_3d(N,N,N);
//		a = minfft_mkaux_t4_3d(N,N,N);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				minfft_dct2(x,y,a);
//				minfft_dst2(x,y,a);
//				minfft_dct3(x,y,a);
//				minfft_dst3(x,y,a);
//				minfft_dct4(x,y,a);
//				minfft_dst4(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		minfft_free_aux(a);
	}
#endif

// performance test of three-dimensional FFTW
#if 0
#include <fftw3.h>
	const int MINT=1;
	const int MAXN=128;
	const int R=10;
	fftw_plan p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		fftw_complex *x = fftw_malloc(N*N*N*sizeof(fftw_complex));
		fftw_complex *y = fftw_malloc(N*N*N*sizeof(fftw_complex));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		// prepare plan
//		p = fftw_plan_dft_3d(N,N,N,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
//		p = fftw_plan_dft_3d(N,N,N,x,y,FFTW_BACKWARD,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
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
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = fftw_malloc(N*N*N*sizeof(double));
		double *y = fftw_malloc(N*N*N*sizeof(double));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = fftw_plan_r2r_3d(N,N,N,x,y,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_3d(N,N,N,x,y,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_3d(N,N,N,x,y,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_3d(N,N,N,x,y,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_3d(N,N,N,x,y,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_3d(N,N,N,x,y,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
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
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		fftw_destroy_plan(p);
	}
#endif

// performance test of three-dimensional Kiss FFT
#if 0
#include "kiss_fftnd.h"
	const int MINT=1;
	const int MAXN=128;
	const int R=10;
	kiss_fftnd_cfg cfg;
	int N,n,Ns[3];
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
		// complex transforms
		kiss_fft_cpx *x = malloc(N*N*N*sizeof(kiss_fft_cpx));
		kiss_fft_cpx *y = malloc(N*N*N*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N*N*N; ++n) {
			x[n].r = (double)rand()/RAND_MAX-0.5;
			x[n].i = (double)rand()/RAND_MAX-0.5;
		}
		// prepare config
		Ns[0] = Ns[1] = Ns[2] = N;
//		cfg = kiss_fftnd_alloc(Ns,3,0,NULL,NULL);
//		cfg = kiss_fftnd_alloc(Ns,3,1,NULL,NULL);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftnd(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,hypot(y[0].r,y[0].i));
		free(x);
		free(y);
		free(cfg);
	}
#endif
}
