#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#include "../dft.h"

int
main (void) {

// one-dimensional DFTs

// compare precision of one-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXBLK=65536*16;
	int N,n;
	double d,dmax; // maximum absolute error
	fftw_plan p; // plan
	struct aux *a; // aux data
	for (N=1; N<=MAXBLK; N*=2) {
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
		a = mkaux_complex(1,&N);
		dft(x,y,a);
		p = fftw_plan_dft_1d(N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		a = mkaux_complex(1,&N);
		idft(x,y,a);
		p = fftw_plan_dft_1d(N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(cabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
#endif
#if 0
		// real DFT
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
		a = mkaux_real_1d(N);
		realdft_1d(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_R2HC,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N/2; ++n) {
			d = log10(fabs(y[2*n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
		if (N>1) {
			d = log10(fabs(y[1]-w[N/2]));
			dmax = (d>dmax)?d:dmax;
		}
		for (n=1; n<N/2; ++n) {
			d = log10(fabs(y[2*n+1]-w[N-n]));
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
		a = mkaux_real_1d(N);
		irealdft_1d(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_HC2R,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		a = mkaux_t2t3(1,&N);
		dct2(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT10,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		a = mkaux_t2t3(1,&N);
		dst2(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT10,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(-y[N-1-n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		a = mkaux_t2t3(1,&N);
		dct3(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT01,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]/2-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
			z[N-1-n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		a = mkaux_t2t3(1,&N);
		dst3(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT01,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(-y[n]/2-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		a = mkaux_t4(1,&N);
		dct4(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_REDFT11,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(y[n]/2-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		a = mkaux_t4(1,&N);
		dst4(x,y,a);
		p = fftw_plan_r2r_1d(N,z,w,FFTW_RODFT11,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N; ++n) {
			d = log10(fabs(-y[n]/2-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(y);
		free(z);
		free(w);
		free_aux(a);
		fftw_destroy_plan(p);
		printf("%8d %g\n",N,dmax);
	}
#endif

// performance test of one-dimensional transforms
#if 0
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	struct aux *a;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
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
		a = mkaux_complex(1,&N);
		// do tests
		T = MINT*MAXBLK/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				dft(x,y,a);
//				idft(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = malloc(N*sizeof(double));
		double *y = malloc(N*sizeof(double));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare aux data
//		a = mkaux_real_1d(N);
//		a = mkaux_t2t3(1,&N);
//		a = mkaux_t4(1,&N);
		// do tests
		T = MAXBLK*MINT/N;
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				realdft_1d(x,y,a);
//				irealdft_1d(x,y,a);
//				dct2(x,y,a);
//				dst2(x,y,a);
//				dct3(x,y,a);
//				dst3(x,y,a);
//				dct4(x,y,a);
//				dst4(x,y,a);
			gettimeofday(&t2,NULL);
			d = (t2.tv_sec-t1.tv_sec)*1000000+(t2.tv_usec-t1.tv_usec);
			v = log2(d/T);
			s += v;
			q += v*v;
		}
		avg = s/R;
		stdd = sqrt((q-s*s/R)/(R-1));
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		free_aux(a);
	}
#endif

// performance test of one-dimensional FFTW
#if 0
#include <fftw3.h>
	const int MINT=10;
	const int MAXBLK=65536*16;
	const int R=100; // repeats
	fftw_plan p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
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
		T = MINT*MAXBLK/N;
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
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,cabs(y[0]));
#endif
#if 0
		// real transforms
		double *x = fftw_malloc(N*sizeof(double));
		double *y = fftw_malloc(N*sizeof(double));
		// prepare test vector
		for (n=0; n<N; ++n)
			x[n] = (double)rand()/RAND_MAX-0.5;
		// prepare plan
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_R2HC,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_HC2R,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT10,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT01,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_REDFT11,FFTW_ESTIMATE);
//		p = fftw_plan_r2r_1d(N,x,y,FFTW_RODFT11,FFTW_ESTIMATE);
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
		printf("%8d %10d %10.3f %10.3f\t\t%g\n",N,R*T,avg,stdd,y[0]);
#endif
		free(x);
		free(y);
		fftw_destroy_plan(p);
	}
#endif

// two-dimensional DFTs

// compare precision of two-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXBLK=1024;
	int N,n,ns[2];
	double d,dmax; // maximum absolute error
	fftw_plan p; // plan
	struct aux *a; // aux data
	for (N=1; N<=MAXBLK; N*=2) {
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_complex(2,ns);
		dft(x,y,a);
		p = fftw_plan_dft_2d(N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(cabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_complex(2,ns);
		idft(x,y,a);
		p = fftw_plan_dft_2d(N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(cabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_t2t3(2,ns);
		dct2(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(fabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_t2t3(2,ns);
		dst2(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(fabs(y[N*N-1-n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_t2t3(2,ns);
		dct3(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(fabs(y[n]/4-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
			z[N*N-1-n] = x[n] = (double)rand()/RAND_MAX-0.5;
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_t2t3(2,ns);
		dst3(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(fabs(y[n]/4-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_t4(2,ns);
		dct4(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(fabs(y[n]/4-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=N;
		// do transforms
		a = mkaux_t4(2,ns);
		dst4(x,y,a);
		p = fftw_plan_r2r_2d(N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N; ++n) {
			d = log10(fabs(y[n]/4-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(y);
		free(z);
		free(w);
		free_aux(a);
		fftw_destroy_plan(p);
		printf("%5d**2 %g\n",N,dmax);
	}
#endif

// performance test of two-dimensional transforms
#if 0
	const int MINT=10;
	const int MAXBLK=1024;
	const int R=100; // repeats
	struct aux *a;
	int N,n,ns[2];
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
#if 0
		// complex transforms
		double complex *x = malloc(N*N*sizeof(double complex));
		double complex *y = malloc(N*N*sizeof(double complex));
		// prepare test vector
		for (n=0; n<N*N; ++n)
			x[n] = \
			(double)rand()/RAND_MAX-0.5+ \
			I*((double)rand()/RAND_MAX-0.5);
		ns[0]=ns[1]=N;
		// prepare aux data
		a = mkaux_complex(2,ns);
		// do tests
		T = MINT*(MAXBLK*MAXBLK)/(N*N);
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				dft(x,y,a);
//				idft(x,y,a);
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
		ns[0]=ns[1]=N;
		// prepare aux data
//		a = mkaux_t2t3(2,ns);
//		a = mkaux_t4(2,ns);
		// do tests
		T = MINT*(MAXBLK*MAXBLK)/(N*N);
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				dct2(x,y,a);
//				dst2(x,y,a);
//				dct3(x,y,a);
//				dst3(x,y,a);
//				dct4(x,y,a);
//				dst4(x,y,a);
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
		free_aux(a);
	}
#endif

// performance test of two-dimensional FFTW
#if 0
#include <fftw3.h>
	const int MINT=10;
	const int MAXBLK=1024;
	const int R=100; // repeats
	fftw_plan p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
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
		T = MINT*(MAXBLK*MAXBLK)/(N*N);
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
		T = MINT*(MAXBLK*MAXBLK)/(N*N);
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

// three-dimensional DFTs

// compare precision of three-dimensional transforms with FFTW
#if 0
#include <unistd.h>
#include <fftw3.h>
	const int MAXBLK=256;
	int N,n,ns[3];
	double d,dmax; // maximum absolute error
	fftw_plan p; // plan
	struct aux *a; // aux data
	for (N=1; N<=MAXBLK; N*=2) {
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_complex(3,ns);
		dft(x,y,a);
		p = fftw_plan_dft_3d(N,N,N,z,w,FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(cabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_complex(3,ns);
		idft(x,y,a);
		p = fftw_plan_dft_3d(N,N,N,z,w,FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(cabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_t2t3(3,ns);
		dct2(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(fabs(y[n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_t2t3(3,ns);
		dst2(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(fabs(-y[N*N*N-1-n]-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_t2t3(3,ns);
		dct3(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(fabs(y[n]/8-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
			z[N*N*N-1-n] = x[n] = (double)rand()/RAND_MAX-0.5;
		// do transforms
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_t2t3(3,ns);
		dst3(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(fabs(-y[n]/8-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_t4(3,ns);
		dct4(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(fabs(y[n]/8-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_t4(3,ns);
		dst4(x,y,a);
		p = fftw_plan_r2r_3d(N,N,N,z,w,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE);
		fftw_execute(p);
		// compare results
		dmax = -HUGE_VAL;
		for (n=0; n<N*N*N; ++n) {
			d = log10(fabs(-y[n]/8-w[n]));
			dmax = (d>dmax)?d:dmax;
		}
#endif
		free(x);
		free(z);
		free(w);
		free_aux(a);
		fftw_destroy_plan(p);
		printf("%5d**3 %g\n",N,dmax);
	}
#endif

// performance test of three-dimensional transforms
#if 0
	const int MINT=10;
	const int MAXBLK=128;
	const int R=100; // repeats
	struct aux *a;
	int N,n,ns[3];
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
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
		ns[0]=ns[1]=ns[2]=N;
		a = mkaux_complex(3,ns);
		// do tests
		T = MINT*(MAXBLK*MAXBLK*MAXBLK)/(N*N*N);
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				dft(x,y,a);
//				idft(x,y,a);
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
		ns[0]=ns[1]=ns[2]=N;
		// prepare aux data
//		a = mkaux_t2t3(3,ns);
//		a = mkaux_t4(3,ns);
		// do tests
		T = MINT*(MAXBLK*MAXBLK*MAXBLK)/(N*N*N);
		s = q = 0.0;
		for (r=0; r<R; ++r) {
			gettimeofday(&t1,NULL);
			for (t=0; t<T; ++t)
//				dct2(x,y,a);
//				dst2(x,y,a);
//				dct3(x,y,a);
//				dst3(x,y,a);
//				dct4(x,y,a);
//				dst4(x,y,a);
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
		free_aux(a);
	}
#endif

// performance test of three-dimensional FFTW
#if 0
#include <fftw3.h>
	const int MINT=10;
	const int MAXBLK=128;
	const int R=100; // repeats
	fftw_plan p;
	int N,n;
	int r,T,t;
	double d,v,s,q,avg,stdd;
	struct timeval t1,t2;
	for (N=1; N<=MAXBLK; N*=2) {
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
		T = MINT*(MAXBLK*MAXBLK*MAXBLK)/(N*N*N);
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
		T = MINT*(MAXBLK*MAXBLK*MAXBLK)/(N*N*N);
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
}
