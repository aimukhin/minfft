#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

#include "../minfft.h"

#define FFTW(OP) fftw_##OP
//#define FFTW(OP) fftwf_##OP
//#define FFTW(OP) fftwl_##OP

int
main (void) {

// one-dimensional DFTs

// compare results of one-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXN=65536*16;
	int N,n,Ncmp;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	FFTW(plan) p;
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
		Ncmp = N;
#if 0
		// complex DFT
		minfft_cmpl *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(minfft_cmpl));
		w = malloc(N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		p = FFTW(plan_dft_1d)(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		FFTW(execute)(p);
#endif
#if 0
		// inverse complex DFT
		minfft_cmpl *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(minfft_cmpl));
		w = malloc(N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_invdft(x,y,a);
		p = FFTW(plan_dft_1d)(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		FFTW(execute)(p);
#endif
#if 0
		// real DFT
		minfft_real *x,*z;
		minfft_cmpl *y,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc((N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(minfft_real));
		w = malloc((N/2+1)*sizeof(minfft_cmpl));
		Ncmp = N/2+1;
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_1d(N);
		minfft_realdft(x,y,a);
		p = FFTW(plan_dft_r2c_1d)(N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x,*z;
		minfft_real *y,*w;
		x = malloc((N/2+1)*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_real));
		z = malloc((N/2+1)*sizeof(minfft_cmpl));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		z[0] = x[0] = (minfft_real)rand()/RAND_MAX-0.5;
		for (n=1; n<N/2; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		z[N/2] = x[N/2] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_1d(N);
		minfft_invrealdft(x,y,a);
		p = FFTW(plan_dft_c2r_1d)(N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-2
		minfft_real *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dct2(x,y,a);
		p = FFTW(plan_r2r_1d)(N,z,w,FFTW_REDFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-2
		minfft_real *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dst2(x,y,a);
		p = FFTW(plan_r2r_1d)(N,z,w,FFTW_RODFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-3
		minfft_real *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dct3(x,y,a);
		p = FFTW(plan_r2r_1d)(N,z,w,FFTW_REDFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-3
		minfft_real *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dst3(x,y,a);
		p = FFTW(plan_r2r_1d)(N,z,w,FFTW_RODFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-4
		minfft_real *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dct4(x,y,a);
		p = FFTW(plan_r2r_1d)(N,z,w,FFTW_REDFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-4
		minfft_real *x,*y,*z,*w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		w = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dst4(x,y,a);
		p = FFTW(plan_r2r_1d)(N,z,w,FFTW_RODFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<Ncmp; ++n) {
			d = fabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%8d %g\n",N,(double)(dmax/vmax));
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		FFTW(destroy_plan)(p);
	}
#endif

// compare results of one-dimensional transforms with Kiss FFT
#if 0
#include <unistd.h>
	const int MAXN=65536*16;
	int N,n;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
#include "kiss_fft.h"
		// complex DFT
		kiss_fft_cfg cfg; // Kiss FFT config
		minfft_cmpl *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(kiss_fft_cpx));
		w = malloc(N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		cfg = kiss_fft_alloc(N,0,NULL,NULL);
		kiss_fft(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fft.h"
		// inverse complex DFT
		kiss_fft_cfg cfg; // Kiss FFT config
		minfft_cmpl *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(kiss_fft_cpx));
		w = malloc(N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_invdft(x,y,a);
		cfg = kiss_fft_alloc(N,1,NULL,NULL);
		kiss_fft(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftr_cfg cfg; // Kiss FFT config
		minfft_real *x;
		minfft_cmpl *y;
		kiss_fft_scalar *z;
		kiss_fft_cpx *w;
		x = malloc(N*sizeof(minfft_real));
		y = malloc((N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(kiss_fft_scalar));
		w = malloc((N/2+1)*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_1d(N);
		minfft_realdft(x,y,a);
		cfg = kiss_fftr_alloc(N,0,NULL,NULL);
		kiss_fftr(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N/2+1; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftr.h"
		// inverse real DFT
		if (N==1)
			continue;
		kiss_fftr_cfg cfg; // Kiss FFT config
		minfft_cmpl *x;
		minfft_real *y;
		kiss_fft_cpx *z;
		kiss_fft_scalar *w;
		x = malloc((N/2+1)*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_real));
		z = malloc((N/2+1)*sizeof(kiss_fft_cpx));
		w = malloc(N*sizeof(kiss_fft_scalar));
		// init inputs
		srand(getpid());
		z[0].r = x[0] = (minfft_real)rand()/RAND_MAX-0.5;
		z[0].i = 0;
		for (n=1; n<N/2; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		z[N/2].r = x[N/2] = (minfft_real)rand()/RAND_MAX-0.5;
		z[N/2].i = 0;
		// do transforms
		a = minfft_mkaux_realdft_1d(N);
		minfft_invrealdft(x,y,a);
		cfg = kiss_fftr_alloc(N,1,NULL,NULL);
		kiss_fftri(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%8d %g\n",N,(double)(dmax/vmax));
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
// with the identity transform
#if 0
#include <unistd.h>
	const int MAXN=65536*16;
	int N,n;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// forward and inverse complex DFT
		minfft_cmpl *x,*y,*z;
		x = malloc(N*sizeof(minfft_cmpl));
		y = malloc(N*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		minfft_invdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(x[n]-z[n]/N);
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// forward and inverse real DFT
		minfft_real *x,*z;
		minfft_cmpl *y;
		x = malloc(N*sizeof(minfft_real));
		y = malloc((N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_1d(N);
		minfft_realdft(x,y,a);
		minfft_invrealdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(x[n]-z[n]/N);
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		minfft_real *x,*y,*z;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dct2(x,y,a);
		minfft_dct3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		minfft_real *x,*y,*z;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_1d(N);
		minfft_dst2(x,y,a);
		minfft_dst3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-4
		minfft_real *x,*y,*z;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dct4(x,y,a);
		minfft_dct4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-4
		minfft_real *x,*y,*z;
		x = malloc(N*sizeof(minfft_real));
		y = malloc(N*sizeof(minfft_real));
		z = malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_1d(N);
		minfft_dst4(x,y,a);
		minfft_dst4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N; ++n) {
			d = fabs(x[n]-z[n]/(2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%8d %g\n",N,(double)(dmax/vmax));
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
		minfft_cmpl *x = malloc(N*sizeof(minfft_cmpl));
		minfft_cmpl *y = malloc(N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
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
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real DFT
		minfft_real *x = malloc(N*sizeof(minfft_real));
		minfft_cmpl *y = malloc(N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare aux data
		a = minfft_mkaux_realdft_1d(N);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				minfft_realdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x = malloc((N/2+1)*sizeof(minfft_cmpl));
		minfft_real *y = malloc(N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N/2+1; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare aux data
		a = minfft_mkaux_realdft_1d(N);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				minfft_invrealdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real transforms
		minfft_real *x = malloc(N*sizeof(minfft_real));
		minfft_real *y = malloc(N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
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
			N,R*T,avg,stdd,(double)fabs(y[0]));
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
	FFTW(plan) p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		minfft_cmpl *x = fftw_malloc(N*sizeof(minfft_cmpl));
		minfft_cmpl *y = fftw_malloc(N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare plan
//		p = FFTW(plan_dft_1d)(N,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
//		p = FFTW(plan_dft_1d)(N,x,y,FFTW_BACKWARD,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real DFT
		minfft_real *x = fftw_malloc(N*sizeof(minfft_real));
		minfft_cmpl *y = fftw_malloc((N/2+1)*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare plan
		p = FFTW(plan_dft_r2c_1d)(N,x,y,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x = fftw_malloc((N/2+1)*sizeof(minfft_cmpl));
		minfft_real *y = fftw_malloc(N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N/2+1; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare plan
		p = FFTW(plan_dft_c2r_1d)(N,x,y,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real transforms
		minfft_real *x = fftw_malloc(N*sizeof(minfft_real));
		minfft_real *y = fftw_malloc(N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = FFTW(plan_r2r_1d)(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_1d)(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_1d)(N,x,y,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_1d)(N,x,y,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_1d)(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_1d)(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
		free(x);
		free(y);
		FFTW(destroy_plan)(p);
	}
#endif

// performance test of one-dimensional Kiss FFT
#if 0
	const int MINT=1;
	const int MAXN=65536*16;
	const int R=10;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
#include "kiss_fft.h"
		// complex transforms
		kiss_fft_cfg cfg;
		kiss_fft_cpx *x = malloc(N*sizeof(kiss_fft_cpx));
		kiss_fft_cpx *y = malloc(N*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N; ++n) {
			x[n].r = (minfft_real)rand()/RAND_MAX-0.5;
			x[n].i = (minfft_real)rand()/RAND_MAX-0.5;
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
			N,R*T,avg,stdd,(double)hypot(y[0].r,y[0].i));
#endif
#if 0
#include "kiss_fftr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftr_cfg cfg;
		kiss_fft_scalar *x = malloc(N*sizeof(kiss_fft_scalar));
		kiss_fft_cpx *y = malloc((N/2+1)*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare config
		cfg = kiss_fftr_alloc(N,0,NULL,NULL);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftr(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)hypot(y[0].r,y[0].i));
#endif
#if 0
#include "kiss_fftr.h"
		// inverse real DFT
		if (N==1)
			continue;
		kiss_fftr_cfg cfg;
		kiss_fft_cpx *x = malloc((N/2+1)*sizeof(kiss_fft_cpx));
		kiss_fft_scalar *y = malloc(N*sizeof(kiss_fft_scalar));
		// prepare test vector
		x[0].r = (minfft_real)rand()/RAND_MAX-0.5;
		x[0].i = 0;
		for (n=1; n<N/2; ++n) {
			x[n].r = (minfft_real)rand()/RAND_MAX-0.5;
			x[n].i = (minfft_real)rand()/RAND_MAX-0.5;
		}
		x[N/2].r = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare config
		cfg = kiss_fftr_alloc(N,1,NULL,NULL);
		// do tests
		T = MINT*MAXN*log2(MAXN+1)/(N*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftri(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
		free(x);
		free(y);
		free(cfg);
	}
#endif

// two-dimensional DFTs

// compare results of two-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXN=1024;
	int N,n,Ncmp;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	FFTW(plan) p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
		Ncmp = N*N;
#if 0
		// complex DFT
		minfft_cmpl *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(minfft_cmpl));
		w = malloc(N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_dft(x,y,a);
		p = FFTW(plan_dft_2d)(N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		FFTW(execute)(p);
#endif
#if 0
		// inverse complex DFT
		minfft_cmpl *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(minfft_cmpl));
		w = malloc(N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_invdft(x,y,a);
		p = FFTW(plan_dft_2d)(N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		FFTW(execute)(p);
#endif
#if 0
		// real DFT
		minfft_real *x,*z;
		minfft_cmpl *y,*w;
		x = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		w = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		Ncmp = N*(N/2+1);
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_2d(N,N);
		minfft_realdft(x,y,a);
		p = FFTW(plan_dft_r2c_2d)(N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x,*z;
		minfft_real *y,*w;
		x = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			y[n] = (minfft_real)rand()/RAND_MAX-0.5;
		a = minfft_mkaux_realdft_2d(N,N);
		minfft_realdft(y,x,a);
		for (n=0; n<N*(N/2+1); ++n)
			z[n] = x[n];
		// do transforms
		minfft_invrealdft(x,y,a);
		p = FFTW(plan_dft_c2r_2d)(N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-2
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dct2(x,y,a);
		p = FFTW(plan_r2r_2d)(N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-2
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dst2(x,y,a);
		p = FFTW(plan_r2r_2d)(N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-3
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dct3(x,y,a);
		p = FFTW(plan_r2r_2d)(N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-3
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dst3(x,y,a);
		p = FFTW(plan_r2r_2d)(N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-4
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dct4(x,y,a);
		p = FFTW(plan_r2r_2d)(N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-4
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		w = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dst4(x,y,a);
		p = FFTW(plan_r2r_2d)(N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<Ncmp; ++n) {
			d = fabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%5d**2 %g\n",N,(double)(dmax/vmax));
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		FFTW(destroy_plan)(p);
	}
#endif

// compare results of two-dimensional transforms with Kiss FFT
#if 0
#include <unistd.h>
	const int MAXN=1024;
	int N,n,Ns[2];
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
#include "kiss_fftnd.h"
		// complex DFT
		kiss_fftnd_cfg cfg; // Kiss FFT config
		minfft_cmpl *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_dft(x,y,a);
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftnd_alloc(Ns,2,0,NULL,NULL);
		kiss_fftnd(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftnd.h"
		// inverse complex DFT
		kiss_fftnd_cfg cfg; // Kiss FFT config
		minfft_cmpl *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_invdft(x,y,a);
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftnd_alloc(Ns,2,1,NULL,NULL);
		kiss_fftnd(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftndr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg; // Kiss FFT config
		minfft_real *x;
		minfft_cmpl *y;
		kiss_fft_scalar *z;
		kiss_fft_cpx *w;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(kiss_fft_scalar));
		w = malloc(N*(N/2+1)*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_2d(N,N);
		minfft_realdft(x,y,a);
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftndr_alloc(Ns,2,0,NULL,NULL);
		kiss_fftndr(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*(N/2+1); ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftndr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg; // Kiss FFT config
		minfft_cmpl *x;
		minfft_real *y;
		kiss_fft_cpx *z;
		kiss_fft_scalar *w;
		x = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*(N/2+1)*sizeof(kiss_fft_cpx));
		w = malloc(N*N*sizeof(kiss_fft_scalar));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			y[n] = (minfft_real)rand()/RAND_MAX-0.5;
		a = minfft_mkaux_realdft_2d(N,N);
		minfft_realdft(y,x,a);
		for (n=0; n<N*(N/2+1); ++n) {
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		minfft_invrealdft(x,y,a);
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftndr_alloc(Ns,2,1,NULL,NULL);
		kiss_fftndri(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%5d**2 %g\n",N,(double)(dmax/vmax));
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
// with the identity transform
#if 0
#include <unistd.h>
	const int MAXN=1024;
	int N,n;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// forward and inverse complex DFT
		minfft_cmpl *x,*y,*z;
		x = malloc(N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_2d(N,N);
		minfft_dft(x,y,a);
		minfft_invdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(x[n]-z[n]/(N*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// forward and inverse real DFT
		minfft_real *x,*z;
		minfft_cmpl *y;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_2d(N,N);
		minfft_realdft(x,y,a);
		minfft_invrealdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(x[n]-z[n]/(N*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		minfft_real *x,*y,*z;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dct2(x,y,a);
		minfft_dct3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		minfft_real *x,*y,*z;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_2d(N,N);
		minfft_dst2(x,y,a);
		minfft_dst3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-4
		minfft_real *x,*y,*z;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dct4(x,y,a);
		minfft_dct4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-4
		minfft_real *x,*y,*z;
		x = malloc(N*N*sizeof(minfft_real));
		y = malloc(N*N*sizeof(minfft_real));
		z = malloc(N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_2d(N,N);
		minfft_dst4(x,y,a);
		minfft_dst4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%5d**2 %g\n",N,(double)(dmax/vmax));
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
		minfft_cmpl *x = malloc(N*N*sizeof(minfft_cmpl));
		minfft_cmpl *y = malloc(N*N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
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
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real DFT
		minfft_real *x = malloc(N*N*sizeof(minfft_real));
		minfft_cmpl *y = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] =(minfft_real)rand()/RAND_MAX-0.5;
		// prepare aux data
		a = minfft_mkaux_realdft_2d(N,N);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				minfft_realdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x = malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		minfft_real *y = malloc(N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*(N/2+1); ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare aux data
		a = minfft_mkaux_realdft_2d(N,N);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				minfft_invrealdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real transforms
		minfft_real *x = malloc(N*N*sizeof(minfft_real));
		minfft_real *y = malloc(N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] =(minfft_real)rand()/RAND_MAX-0.5;
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
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
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
	FFTW(plan) p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		minfft_cmpl *x = fftw_malloc(N*N*sizeof(minfft_cmpl));
		minfft_cmpl *y = fftw_malloc(N*N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare plan
//		p = FFTW(plan_dft_2d)(N,N,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
//		p = FFTW(plan_dft_2d)(N,N,x,y,FFTW_BACKWARD,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real DFT
		minfft_real *x = fftw_malloc(N*N*sizeof(minfft_real));
		minfft_cmpl *y = fftw_malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare plan
		p = FFTW(plan_dft_r2c_2d)(N,N,x,y,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x = fftw_malloc(N*(N/2+1)*sizeof(minfft_cmpl));
		minfft_real *y = fftw_malloc(N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*(N/2+1); ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare plan
		p = FFTW(plan_dft_c2r_2d)(N,N,x,y,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real transforms
		minfft_real *x = fftw_malloc(N*N*sizeof(minfft_real));
		minfft_real *y = fftw_malloc(N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = FFTW(plan_r2r_2d)(N,N,x,y,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_2d)(N,N,x,y,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_2d)(N,N,x,y,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_2d)(N,N,x,y,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_2d)(N,N,x,y,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_2d)(N,N,x,y,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
		free(x);
		free(y);
		FFTW(destroy_plan)(p);
	}
#endif

// performance test of two-dimensional Kiss FFT
#if 0
	const int MINT=1;
	const int MAXN=1024;
	const int R=10;
	int N,n,Ns[2];
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
#include "kiss_fftnd.h"
		// complex transforms
		kiss_fftnd_cfg cfg;
		kiss_fft_cpx *x = malloc(N*N*sizeof(kiss_fft_cpx));
		kiss_fft_cpx *y = malloc(N*N*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N*N; ++n) {
			x[n].r = (minfft_real)rand()/RAND_MAX-0.5;
			x[n].i = (minfft_real)rand()/RAND_MAX-0.5;
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
			N,R*T,avg,stdd,(double)hypot(y[0].r,y[0].i));
#endif
#if 0
#include "kiss_fftndr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg;
		kiss_fft_scalar *x = malloc(N*N*sizeof(kiss_fft_scalar));
		kiss_fft_cpx *y = malloc(N*(N/2+1)*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare config
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftndr_alloc(Ns,2,0,NULL,NULL);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftndr(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)hypot(y[0].r,y[0].i));
#endif
#if 0
#include "kiss_fftndr.h"
		// inverse real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg;
		kiss_fft_cpx *x = malloc(N*(N/2+1)*sizeof(kiss_fft_cpx));
		kiss_fft_scalar *y = malloc(N*N*sizeof(kiss_fft_scalar));
		// prepare test vector
		for (n=0; n<N*(N/2+1); ++n) {
			x[n].r = (minfft_real)rand()/RAND_MAX-0.5;
			x[n].i = (minfft_real)rand()/RAND_MAX-0.5;
		}
		// prepare config
		Ns[0] = Ns[1] = N;
		cfg = kiss_fftndr_alloc(Ns,2,1,NULL,NULL);
		// do tests
		T = MINT*MAXN*MAXN*2*log2(MAXN+1)/(N*N*2*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftndri(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**2 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
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
	int N,n,Ncmp;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	FFTW(plan) p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
		Ncmp = N*N*N;
#if 0
		// complex DFT
		minfft_cmpl *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(minfft_cmpl));
		w = malloc(N*N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_dft(x,y,a);
		p = FFTW(plan_dft_3d)(N,N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		FFTW(execute)(p);
#endif
#if 0
		// inverse complex DFT
		minfft_cmpl *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(minfft_cmpl));
		w = malloc(N*N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_invdft(x,y,a);
		p = FFTW(plan_dft_3d)(N,N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		FFTW(execute)(p);
#endif
#if 0
		// real DFT
		minfft_real *x,*z;
		minfft_cmpl *y,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		w = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		Ncmp = N*N*(N/2+1);
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_3d(N,N,N);
		minfft_realdft(x,y,a);
		p = FFTW(plan_dft_r2c_3d)(N,N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x,*z;
		minfft_real *y,*w;
		x = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			y[n] = (minfft_real)rand()/RAND_MAX-0.5;
		a = minfft_mkaux_realdft_3d(N,N,N);
		minfft_realdft(y,x,a);
		for (n=0; n<N*N*(N/2+1); ++n)
			z[n] = x[n];
		// do transforms
		minfft_invrealdft(x,y,a);
		p = FFTW(plan_dft_c2r_3d)(N,N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-2
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dct2(x,y,a);
		p = FFTW(plan_r2r_3d)(N,N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-2
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dst2(x,y,a);
		p = FFTW(plan_r2r_3d)(N,N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-3
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dct3(x,y,a);
		p = FFTW(plan_r2r_3d)(N,N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-3
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dst3(x,y,a);
		p = FFTW(plan_r2r_3d)(N,N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DCT-4
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dct4(x,y,a);
		p = FFTW(plan_r2r_3d)(N,N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if 0
		// DST-4
		minfft_real *x,*y,*z,*w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		w = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dst4(x,y,a);
		p = FFTW(plan_r2r_3d)(N,N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<Ncmp; ++n) {
			d = fabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
		printf("%5d**3 %g\n",N,(double)(dmax/vmax));
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
	const int MAXN=128;
	int N,n,Ns[3];
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
#include "kiss_fftnd.h"
		// complex DFT
		kiss_fftnd_cfg cfg; // Kiss FFT config
		minfft_cmpl *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_dft(x,y,a);
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftnd_alloc(Ns,3,0,NULL,NULL);
		kiss_fftnd(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftnd.h"
		// inverse complex DFT
		kiss_fftnd_cfg cfg; // Kiss FFT config
		minfft_cmpl *x,*y;
		kiss_fft_cpx *z,*w;
		x = malloc(N*N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(kiss_fft_cpx));
		w = malloc(N*N*N*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n) {
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_invdft(x,y,a);
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftnd_alloc(Ns,3,1,NULL,NULL);
		kiss_fftnd(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftndr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg; // Kiss FFT config
		minfft_real *x;
		minfft_cmpl *y;
		kiss_fft_scalar *z;
		kiss_fft_cpx *w;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(kiss_fft_scalar));
		w = malloc(N*N*(N/2+1)*sizeof(kiss_fft_cpx));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			z[n] = x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_3d(N,N,N);
		minfft_realdft(x,y,a);
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftndr_alloc(Ns,3,0,NULL,NULL);
		kiss_fftndr(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*(N/2+1); ++n) {
			d = fabs(y[n]-(w[n].r+I*w[n].i));
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
#include "kiss_fftndr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg; // Kiss FFT config
		minfft_cmpl *x;
		minfft_real *y;
		kiss_fft_cpx *z;
		kiss_fft_scalar *w;
		x = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*(N/2+1)*sizeof(kiss_fft_cpx));
		w = malloc(N*N*N*sizeof(kiss_fft_scalar));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			y[n] = (minfft_real)rand()/RAND_MAX-0.5;
		a = minfft_mkaux_realdft_3d(N,N,N);
		minfft_realdft(y,x,a);
		for (n=0; n<N*N*(N/2+1); ++n) {
			z[n].r = creal(x[n]);
			z[n].i = cimag(x[n]);
		}	
		// do transforms
		minfft_invrealdft(x,y,a);
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftndr_alloc(Ns,3,1,NULL,NULL);
		kiss_fftndri(cfg,z,w);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(y[n]-w[n]);
			dmax = (d>dmax)?d:dmax;
			v = fabs(y[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%5d**3 %g\n",N,(double)(dmax/vmax));
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
// with the identity transform
#if 0
#include <unistd.h>
	const int MAXN=128;
	int N,n;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// forward and inverse complex DFT
		minfft_cmpl *x,*y,*z;
		x = malloc(N*N*N*sizeof(minfft_cmpl));
		y = malloc(N*N*N*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
		a = minfft_mkaux_dft_3d(N,N,N);
		minfft_dft(x,y,a);
		minfft_invdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(x[n]-z[n]/(N*N*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// forward and inverse real DFT
		minfft_real *x,*z;
		minfft_cmpl *y;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		z = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_realdft_3d(N,N,N);
		minfft_realdft(x,y,a);
		minfft_invrealdft(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(x[n]-z[n]/(N*N*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-2 and DCT-3
		minfft_real *x,*y,*z;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dct2(x,y,a);
		minfft_dct3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-2 and DST-3
		minfft_real *x,*y,*z;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t2t3_3d(N,N,N);
		minfft_dst2(x,y,a);
		minfft_dst3(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DCT-4
		minfft_real *x,*y,*z;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dct4(x,y,a);
		minfft_dct4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
#if 0
		// DST-4
		minfft_real *x,*y,*z;
		x = malloc(N*N*N*sizeof(minfft_real));
		y = malloc(N*N*N*sizeof(minfft_real));
		z = malloc(N*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a = minfft_mkaux_t4_3d(N,N,N);
		minfft_dst4(x,y,a);
		minfft_dst4(y,z,a);
		// compare results
		dmax = 0;
		vmax = 0;
		for (n=0; n<N*N*N; ++n) {
			d = fabs(x[n]-z[n]/(2*N*2*N*2*N));
			dmax = (d>dmax)?d:dmax;
			v = fabs(x[n]);
			vmax = (v>vmax)?v:vmax;
		}
#endif
		// print results
		printf("%5d**3 %g\n",N,(double)(dmax/vmax));
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
		minfft_cmpl *x = malloc(N*N*N*sizeof(minfft_cmpl));
		minfft_cmpl *y = malloc(N*N*N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
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
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real DFT
		minfft_real *x = malloc(N*N*N*sizeof(minfft_real));
		minfft_cmpl *y = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare aux data
		a = minfft_mkaux_realdft_3d(N,N,N);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				minfft_realdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x = malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		minfft_real *y = malloc(N*N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*N*(N/2+1); ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare aux data
		a = minfft_mkaux_realdft_3d(N,N,N);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				minfft_invrealdft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real transforms
		minfft_real *x = malloc(N*N*N*sizeof(minfft_real));
		minfft_real *y = malloc(N*N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] =(minfft_real)rand()/RAND_MAX-0.5;
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
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
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
	FFTW(plan) p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
		// complex transforms
		minfft_cmpl *x = fftw_malloc(N*N*N*sizeof(minfft_cmpl));
		minfft_cmpl *y = fftw_malloc(N*N*N*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare plan
//		p = FFTW(plan_dft_3d)(N,N,N,x,y,FFTW_FORWARD,FFTW_ESTIMATE);
//		p = FFTW(plan_dft_3d)(N,N,N,x,y,FFTW_BACKWARD,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real DFT
		minfft_real *x = fftw_malloc(N*N*N*sizeof(minfft_real));
		minfft_cmpl *y = fftw_malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare plan
		p = FFTW(plan_dft_r2c_3d)(N,N,N,x,y,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// inverse real DFT
		minfft_cmpl *x = fftw_malloc(N*N*(N/2+1)*sizeof(minfft_cmpl));
		minfft_real *y = fftw_malloc(N*N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*N*(N/2+1); ++n)
			x[n] = \
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// prepare plan
		p = FFTW(plan_dft_c2r_3d)(N,N,N,x,y,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
#if 0
		// real transforms
		minfft_real *x = fftw_malloc(N*N*N*sizeof(minfft_real));
		minfft_real *y = fftw_malloc(N*N*N*sizeof(minfft_real));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = FFTW(plan_r2r_3d)(N,N,N,x,y,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_3d)(N,N,N,x,y,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_3d)(N,N,N,x,y,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_3d)(N,N,N,x,y,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_3d)(N,N,N,x,y,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = FFTW(plan_r2r_3d)(N,N,N,x,y,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				FFTW(execute)(p);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
		free(x);
		free(y);
		FFTW(destroy_plan)(p);
	}
#endif

// performance test of three-dimensional Kiss FFT
#if 0
	const int MINT=1;
	const int MAXN=128;
	const int R=10;
	int N,n,Ns[3];
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXN; N*=2) {
#if 0
#include "kiss_fftnd.h"
		// complex transforms
		kiss_fftnd_cfg cfg;
		kiss_fft_cpx *x = malloc(N*N*N*sizeof(kiss_fft_cpx));
		kiss_fft_cpx *y = malloc(N*N*N*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N*N*N; ++n) {
			x[n].r = (minfft_real)rand()/RAND_MAX-0.5;
			x[n].i = (minfft_real)rand()/RAND_MAX-0.5;
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
			N,R*T,avg,stdd,(double)hypot(y[0].r,y[0].i));
#endif
#if 0
#include "kiss_fftndr.h"
		// real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg;
		kiss_fft_scalar *x = malloc(N*N*N*sizeof(kiss_fft_scalar));
		kiss_fft_cpx *y = malloc(N*N*(N/2+1)*sizeof(kiss_fft_cpx));
		// prepare test vector
		for (n=0; n<N*N*N; ++n)
			x[n] = (minfft_real)rand()/RAND_MAX-0.5;
		// prepare config
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftndr_alloc(Ns,3,0,NULL,NULL);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftndr(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)hypot(y[0].r,y[0].i));
#endif
#if 0
#include "kiss_fftndr.h"
		// inverse real DFT
		if (N==1)
			continue;
		kiss_fftndr_cfg cfg;
		kiss_fft_cpx *x = malloc(N*N*(N/2+1)*sizeof(kiss_fft_cpx));
		kiss_fft_scalar *y = malloc(N*N*N*sizeof(kiss_fft_scalar));
		// prepare test vector
		for (n=0; n<N*N*(N/2+1); ++n) {
			x[n].r = (minfft_real)rand()/RAND_MAX-0.5;
			x[n].i = (minfft_real)rand()/RAND_MAX-0.5;
		}
		// prepare config
		Ns[0] = Ns[1] = Ns[2] = N;
		cfg = kiss_fftndr_alloc(Ns,3,1,NULL,NULL);
		// do tests
		T = MINT*MAXN*MAXN*MAXN*3*log2(MAXN+1)/(N*N*N*3*log2(N+1));
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
				kiss_fftndri(cfg,x,y);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%5d**3 %10d %10.3f %10.3f\t\t%g\n",
			N,R*T,avg,stdd,(double)fabs(y[0]));
#endif
		free(x);
		free(y);
		free(cfg);
	}
#endif
}
