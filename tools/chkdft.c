// Compare minfft with FFTW
// Dimensionality: D1 D2 D3
// Transforms: DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4

// macros for calling variants of FFTW
#ifndef FFTW_SFX
#define FFTW_SFX
#endif
#define FFTW(F) FFTW2(F,FFTW_SFX)
#define FFTW2(F,S) FFTW22(F,S)
#define FFTW22(F,S) fftw##S##_##F

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fftw3.h>
#include <math.h>
#include <tgmath.h>

#include "minfft.h"

int
main (void) {

// one-dimensional DFTs
#if D1
	const int MAXN=65536*16;
	int N,n,Ncmp;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_real r; // relative error
	FFTW(plan) p;
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
		Ncmp=N;
#if DFT || INVDFT
		// complex transforms
		minfft_cmpl *x,*y,*z,*w;
		x=malloc(N*sizeof(minfft_cmpl));
		y=malloc(N*sizeof(minfft_cmpl));
		z=malloc(N*sizeof(minfft_cmpl));
		w=malloc(N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n]=x[n]=\
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
#if DFT
		a=minfft_mkaux_dft_1d(N);
		minfft_dft(x,y,a);
		p=FFTW(plan_dft_1d)(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif INVDFT
		a=minfft_mkaux_dft_1d(N);
		minfft_invdft(x,y,a);
		p=FFTW(plan_dft_1d)(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#endif
#if REALDFT
		// real DFT
		minfft_real *x,*z;
		minfft_cmpl *y,*w;
		x=malloc(N*sizeof(minfft_real));
		y=malloc((N/2+1)*sizeof(minfft_cmpl));
		z=malloc(N*sizeof(minfft_real));
		w=malloc((N/2+1)*sizeof(minfft_cmpl));
		Ncmp=N/2+1;
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n]=x[n]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a=minfft_mkaux_realdft_1d(N);
		minfft_realdft(x,y,a);
		p=FFTW(plan_dft_r2c_1d)(N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if INVREALDFT
		// inverse real DFT
		minfft_cmpl *x,*z;
		minfft_real *y,*w;
		x=malloc((N/2+1)*sizeof(minfft_cmpl));
		y=malloc(N*sizeof(minfft_real));
		z=malloc((N/2+1)*sizeof(minfft_cmpl));
		w=malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		z[0]=x[0]=(minfft_real)rand()/RAND_MAX-0.5;
		for (n=1; n<N/2; ++n)
			z[n]=x[n]=\
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		z[N/2]=x[N/2]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a=minfft_mkaux_realdft_1d(N);
		minfft_invrealdft(x,y,a);
		p=FFTW(plan_dft_c2r_1d)(N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
		// real symmetric transforms
		minfft_real *x,*y,*z,*w;
		x=malloc(N*sizeof(minfft_real));
		y=malloc(N*sizeof(minfft_real));
		z=malloc(N*sizeof(minfft_real));
		w=malloc(N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<N; ++n)
			z[n]=x[n]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
#if DCT2
		a=minfft_mkaux_t2t3_1d(N);
		minfft_dct2(x,y,a);
		p=FFTW(plan_r2r_1d)(N,z,w,FFTW_REDFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST2
		a=minfft_mkaux_t2t3_1d(N);
		minfft_dst2(x,y,a);
		p=FFTW(plan_r2r_1d)(N,z,w,FFTW_RODFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DCT3
		a=minfft_mkaux_t2t3_1d(N);
		minfft_dct3(x,y,a);
		p=FFTW(plan_r2r_1d)(N,z,w,FFTW_REDFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST3
		a=minfft_mkaux_t2t3_1d(N);
		minfft_dst3(x,y,a);
		p=FFTW(plan_r2r_1d)(N,z,w,FFTW_RODFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DCT4
		a=minfft_mkaux_t4_1d(N);
		minfft_dct4(x,y,a);
		p=FFTW(plan_r2r_1d)(N,z,w,FFTW_REDFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST4
		a=minfft_mkaux_t4_1d(N);
		minfft_dst4(x,y,a);
		p=FFTW(plan_r2r_1d)(N,z,w,FFTW_RODFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#endif
		// compare results
		dmax=0;
		vmax=0;
		for (n=0; n<Ncmp; ++n) {
			d=fabs(y[n]-w[n]);
			dmax=(d>dmax)?d:dmax;
			v=fabs(y[n]);
			vmax=(v>vmax)?v:vmax;
		}
		r=dmax/vmax;
		printf("%12d\t%g\n",N,(double)r);
		if (r>=T)
			// threshold exceeded
			return 1;
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		FFTW(destroy_plan)(p);
	}
#endif

// two-dimensional DFTs
#if D2
	const int MAXN=1024;
	int N,n,Ncmp;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_real r; // relative error
	FFTW(plan) p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
		Ncmp=2*N*N;
#if DFT || INVDFT
		// complex transforms
		minfft_cmpl *x,*y,*z,*w;
		x=malloc(2*N*N*sizeof(minfft_cmpl));
		y=malloc(2*N*N*sizeof(minfft_cmpl));
		z=malloc(2*N*N*sizeof(minfft_cmpl));
		w=malloc(2*N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<2*N*N; ++n)
			z[n]=x[n]=\
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
#if DFT
		// do transforms
		a=minfft_mkaux_dft_2d(2*N,N);
		minfft_dft(x,y,a);
		p=FFTW(plan_dft_2d)(2*N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif INVDFT
		a=minfft_mkaux_dft_2d(2*N,N);
		minfft_invdft(x,y,a);
		p=FFTW(plan_dft_2d)(2*N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#endif
#if REALDFT
		// real DFT
		minfft_real *x,*z;
		minfft_cmpl *y,*w;
		x=malloc(2*N*N*sizeof(minfft_real));
		z=malloc(2*N*N*sizeof(minfft_real));
		y=malloc(2*N*(N/2+1)*sizeof(minfft_cmpl));
		w=malloc(2*N*(N/2+1)*sizeof(minfft_cmpl));
		Ncmp=2*N*(N/2+1);
		// init inputs
		srand(getpid());
		for (n=0; n<2*N*N; ++n)
			z[n]=x[n]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a=minfft_mkaux_realdft_2d(2*N,N);
		minfft_realdft(x,y,a);
		p=FFTW(plan_dft_r2c_2d)(2*N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if INVREALDFT
		// inverse real DFT
		minfft_cmpl *x,*z;
		minfft_real *y,*w;
		x=malloc(2*N*(N/2+1)*sizeof(minfft_cmpl));
		z=malloc(2*N*(N/2+1)*sizeof(minfft_cmpl));
		y=malloc(2*N*N*sizeof(minfft_real));
		w=malloc(2*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<2*N*N; ++n)
			y[n]=(minfft_real)rand()/RAND_MAX-0.5;
		a=minfft_mkaux_realdft_2d(2*N,N);
		minfft_realdft(y,x,a);
		for (n=0; n<2*N*(N/2+1); ++n)
			z[n]=x[n];
		// do transforms
		minfft_invrealdft(x,y,a);
		p=FFTW(plan_dft_c2r_2d)(2*N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
		// real symmetric transforms
		minfft_real *x,*y,*z,*w;
		x=malloc(2*N*N*sizeof(minfft_real));
		y=malloc(2*N*N*sizeof(minfft_real));
		z=malloc(2*N*N*sizeof(minfft_real));
		w=malloc(2*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<2*N*N; ++n)
			z[n]=x[n]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
#if DCT2
		a=minfft_mkaux_t2t3_2d(2*N,N);
		minfft_dct2(x,y,a);
		p=FFTW(plan_r2r_2d)(2*N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST2
		a=minfft_mkaux_t2t3_2d(2*N,N);
		minfft_dst2(x,y,a);
		p=FFTW(plan_r2r_2d)(2*N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DCT3
		a=minfft_mkaux_t2t3_2d(2*N,N);
		minfft_dct3(x,y,a);
		p=FFTW(plan_r2r_2d)(2*N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST3
		a=minfft_mkaux_t2t3_2d(2*N,N);
		minfft_dst3(x,y,a);
		p=FFTW(plan_r2r_2d)(2*N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DCT4
		a=minfft_mkaux_t4_2d(2*N,N);
		minfft_dct4(x,y,a);
		p=FFTW(plan_r2r_2d)(2*N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST4
		a=minfft_mkaux_t4_2d(2*N,N);
		minfft_dst4(x,y,a);
		p=FFTW(plan_r2r_2d)(2*N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#endif
		// compare results
		dmax=0;
		vmax=0;
		for (n=0; n<Ncmp; ++n) {
			d=fabs(y[n]-w[n]);
			dmax=(d>dmax)?d:dmax;
			v=fabs(y[n]);
			vmax=(v>vmax)?v:vmax;
		}
		// print results
		char s[256];
		snprintf(s,256,"%d*%d",2*N,N);
		r=dmax/vmax;
		printf("%12s\t%g\n",s,(double)r);
		if (r>=T)
			// threshold exceeded
			return 1;
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
		FFTW(destroy_plan)(p);
	}
#endif

// three-dimensional DFTs
#if D3
	const int MAXN=64;
	int N,n,Ncmp;
	minfft_real d,dmax,v,vmax; // current and maximum absolute values
	minfft_real r; // relative error
	FFTW(plan) p; // plan
	minfft_aux *a; // aux data
	for (N=1; N<=MAXN; N*=2) {
		Ncmp=4*N*2*N*N;
#if DFT || INVDFT
		// complex transforms
		minfft_cmpl *x,*y,*z,*w;
		x=malloc(4*N*2*N*N*sizeof(minfft_cmpl));
		y=malloc(4*N*2*N*N*sizeof(minfft_cmpl));
		z=malloc(4*N*2*N*N*sizeof(minfft_cmpl));
		w=malloc(4*N*2*N*N*sizeof(minfft_cmpl));
		// init inputs
		srand(getpid());
		for (n=0; n<4*N*2*N*N; ++n)
			z[n]=x[n]=\
			(minfft_real)rand()/RAND_MAX-0.5+ \
			I*((minfft_real)rand()/RAND_MAX-0.5);
		// do transforms
#if DFT
		a=minfft_mkaux_dft_3d(4*N,2*N,N);
		minfft_dft(x,y,a);
		p=FFTW(plan_dft_3d)(4*N,2*N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif INVDFT
		a=minfft_mkaux_dft_3d(4*N,2*N,N);
		minfft_invdft(x,y,a);
		p=FFTW(plan_dft_3d)(4*N,2*N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#endif
#if REALDFT
		// real DFT
		minfft_real *x,*z;
		minfft_cmpl *y,*w;
		x=malloc(4*N*2*N*N*sizeof(minfft_real));
		z=malloc(4*N*2*N*N*sizeof(minfft_real));
		y=malloc(4*N*2*N*(N/2+1)*sizeof(minfft_cmpl));
		w=malloc(4*N*2*N*(N/2+1)*sizeof(minfft_cmpl));
		Ncmp=4*N*2*N*(N/2+1);
		// init inputs
		srand(getpid());
		for (n=0; n<4*N*2*N*N; ++n)
			z[n]=x[n]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
		a=minfft_mkaux_realdft_3d(4*N,2*N,N);
		minfft_realdft(x,y,a);
		p=FFTW(plan_dft_r2c_3d)(4*N,2*N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if INVREALDFT
		// inverse real DFT
		minfft_cmpl *x,*z;
		minfft_real *y,*w;
		x=malloc(4*N*2*N*(N/2+1)*sizeof(minfft_cmpl));
		z=malloc(4*N*2*N*(N/2+1)*sizeof(minfft_cmpl));
		y=malloc(4*N*2*N*N*sizeof(minfft_real));
		w=malloc(4*N*2*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<4*N*2*N*N; ++n)
			y[n]=(minfft_real)rand()/RAND_MAX-0.5;
		a=minfft_mkaux_realdft_3d(4*N,2*N,N);
		minfft_realdft(y,x,a);
		for (n=0; n<4*N*2*N*(N/2+1); ++n)
			z[n]=x[n];
		// do transforms
		minfft_invrealdft(x,y,a);
		p=FFTW(plan_dft_c2r_3d)(4*N,2*N,N,z,w,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
		// real symmetric transforms
		minfft_real *x,*y,*z,*w;
		x=malloc(4*N*2*N*N*sizeof(minfft_real));
		y=malloc(4*N*2*N*N*sizeof(minfft_real));
		z=malloc(4*N*2*N*N*sizeof(minfft_real));
		w=malloc(4*N*2*N*N*sizeof(minfft_real));
		// init inputs
		srand(getpid());
		for (n=0; n<4*N*2*N*N; ++n)
			z[n]=x[n]=(minfft_real)rand()/RAND_MAX-0.5;
		// do transforms
#if DCT2
		a=minfft_mkaux_t2t3_3d(4*N,2*N,N);
		minfft_dct2(x,y,a);
		p=FFTW(plan_r2r_3d)(4*N,2*N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST2
		a=minfft_mkaux_t2t3_3d(4*N,2*N,N);
		minfft_dst2(x,y,a);
		p=FFTW(plan_r2r_3d)(4*N,2*N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DCT3
		a=minfft_mkaux_t2t3_3d(4*N,2*N,N);
		minfft_dct3(x,y,a);
		p=FFTW(plan_r2r_3d)(4*N,2*N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST3
		a=minfft_mkaux_t2t3_3d(4*N,2*N,N);
		minfft_dst3(x,y,a);
		p=FFTW(plan_r2r_3d)(4*N,2*N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DCT4
		a=minfft_mkaux_t4_3d(4*N,2*N,N);
		minfft_dct4(x,y,a);
		p=FFTW(plan_r2r_3d)(4*N,2*N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#elif DST4
		a=minfft_mkaux_t4_3d(4*N,2*N,N);
		minfft_dst4(x,y,a);
		p=FFTW(plan_r2r_3d)(4*N,2*N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		FFTW(execute)(p);
#endif
#endif
		// compare results
		dmax=0;
		vmax=0;
		for (n=0; n<Ncmp; ++n) {
			d=fabs(y[n]-w[n]);
			dmax=(d>dmax)?d:dmax;
			v=fabs(y[n]);
			vmax=(v>vmax)?v:vmax;
		}
		// print results
		char s[256];
		snprintf(s,256,"%d*%d*%d",4*N,2*N,N);
		r=dmax/vmax;
		printf("%12s\t%g\n",s,(double)r);
		if (r>=T)
			// threshold exceeded
			return 1;
		// free resources
		free(x);
		free(y);
		free(z);
		free(w);
		minfft_free_aux(a);
	}
#endif

return 0;

}
