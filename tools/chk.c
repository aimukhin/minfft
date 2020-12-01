#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <tgmath.h>

#define D_MAX 4
#define T_MAX (1<<30)
#define US_MIN 100000

#include "minfft.h"

#include <fftw3.h>
#ifndef FFTW_SFX
#define FFTW_SFX
#endif
#define FFTW(F) FFTW2(F,FFTW_SFX)
#define FFTW2(F,S) FFTW22(F,S)
#define FFTW22(F,S) fftw##S##_##F

struct {
	int d,n[D_MAX];
	void *in1,*in2;
	void *out1,*out2;
	void *x1,*x2;
	int insz,outsz;
} P;

int
parse (char *s) {
	char *h,*t;
	int v;
	h=t=s;
	P.d=0;
	while (1) {
		if (P.d==D_MAX)
			return -1;
		while (isdigit(*t)) ++t;
		if (*t=='x') {
			*t='\0';
			v=atoi(h);
			if (v&(v-1))
				return -2;
			++P.d;
			P.n[P.d-1]=v;
			++t;
			h=t;
		} else if (*t=='\0') {
			v=atoi(h);
			if (v&(v-1))
				return -2;
			++P.d;
			P.n[P.d-1]=v;
			return 0;
		} else
			return -3;
	}
}

void
alloc (void) {
	int i,p;
#if DFT||INVDFT
	for (p=1,i=0; i<P.d; p*=P.n[i++]);
	P.insz=2*p;
	P.outsz=2*p;
#elif REALDFT
	for (p=1,i=0; i<P.d-1; p*=P.n[i++]);
	P.insz=p*P.n[P.d-1];
	P.outsz=p*(P.n[P.d-1]+2);
#elif INVREALDFT
	for (p=1,i=0; i<P.d-1; p*=P.n[i++]);
	P.insz=p*(P.n[P.d-1]+2);
	P.outsz=p*P.n[P.d-1];
#elif DCT2||DST2||DCT3||DST3||DCT4||DST4
	for (p=1,i=0; i<P.d; p*=P.n[i++]);
	P.insz=p;
	P.outsz=p;
#endif
	P.in1=malloc(P.insz*sizeof(minfft_real));
	P.in2=malloc(P.insz*sizeof(minfft_real));
	P.out1=malloc(P.outsz*sizeof(minfft_real));
	P.out2=malloc(P.outsz*sizeof(minfft_real));
}

void
init (void) {
	int i;
	minfft_real *in1=P.in1,*in2=P.in2;
#if !INVREALDFT
	for (i=0; i<P.insz; ++i)
		in1[i]=in2[i]=(minfft_real)rand()/RAND_MAX-0.5;
#else
	// tricky
	minfft_aux *a=minfft_mkaux_realdft(P.d,P.n);
	minfft_real *x=P.out1;
	for (i=0; i<P.outsz; ++i)
		x[i]=(minfft_real)rand()/RAND_MAX-0.5;
	minfft_realdft(x,P.in1,a);
	minfft_free_aux(a);
	for (i=0; i<P.insz; ++i)
		in2[i]=in1[i];
#endif
}

void
setup_minfft (void) {
#if DFT||INVDFT
	P.x1=minfft_mkaux_dft(P.d,P.n);
#elif REALDFT||INVREALDFT
	P.x1=minfft_mkaux_realdft(P.d,P.n);
#elif DCT2||DST2||DCT3||DST3
	P.x1=minfft_mkaux_t2t3(P.d,P.n);
#elif DCT4||DST4
	P.x1=minfft_mkaux_t4(P.d,P.n);
#endif
}

void
doit_minfft (int T) {
	int t;
	for (t=0; t<T; ++t)
#if DFT
		minfft_dft(P.in1,P.out1,P.x1);
#elif INVDFT
		minfft_invdft(P.in1,P.out1,P.x1);
#elif REALDFT
		minfft_realdft(P.in1,P.out1,P.x1);
#elif INVREALDFT
		minfft_invrealdft(P.in1,P.out1,P.x1);
#elif DCT2
		minfft_dct2(P.in1,P.out1,P.x1);
#elif DST2
		minfft_dst2(P.in1,P.out1,P.x1);
#elif DCT3
		minfft_dct3(P.in1,P.out1,P.x1);
#elif DST3
		minfft_dst3(P.in1,P.out1,P.x1);
#elif DCT4
		minfft_dct4(P.in1,P.out1,P.x1);
#elif DST4
		minfft_dst4(P.in1,P.out1,P.x1);
#endif
}

void
cleanup_minfft (void) {
	minfft_free_aux(P.x1);
}

void
setup_fftw (void) {
#if DFT
	P.x2=FFTW(plan_dft)(P.d,P.n,P.in2,P.out2,FFTW_FORWARD,FFTW_ESTIMATE);
#elif INVDFT
	P.x2=FFTW(plan_dft)(P.d,P.n,P.in2,P.out2,FFTW_BACKWARD,FFTW_ESTIMATE);
#elif REALDFT
	P.x2=FFTW(plan_dft_r2c)(P.d,P.n,P.in2,P.out2,FFTW_ESTIMATE);
#elif INVREALDFT
	P.x2=FFTW(plan_dft_c2r)(P.d,P.n,P.in2,P.out2,FFTW_ESTIMATE);
#else
	fftw_r2r_kind kinds[D_MAX];
	for (int i=0; i<P.d; ++i)
#if DCT2
		kinds[i]=FFTW_REDFT10;
#elif DST2
		kinds[i]=FFTW_RODFT10;
#elif DCT3
		kinds[i]=FFTW_REDFT01;
#elif DST3
		kinds[i]=FFTW_RODFT01;
#elif DCT4
		kinds[i]=FFTW_REDFT11;
#elif DST4
		kinds[i]=FFTW_RODFT11;
#endif
	P.x2=FFTW(plan_r2r)(P.d,P.n,P.in2,P.out2,kinds,FFTW_ESTIMATE);
#endif
}

void
doit_fftw (int T) {
	int t;
	for (t=0; t<T; ++t)
		FFTW(execute)(P.x2);
}

void
cleanup_fftw (void) {
	FFTW(destroy_plan)(P.x2);
}

void
dealloc (void) {
	free(P.in1);
	free(P.in2);
	free(P.out1);
	free(P.out2);
}

double
cmp (void) {
	double d,v,dmax,vmax;
	minfft_real *o1,*o2;
	int i;
	doit_minfft(1);
	doit_fftw(1);
	dmax=0;
	vmax=0;
	o1=P.out1;
	o2=P.out2;
	for (i=0; i<P.outsz; ++i) {
		d=fabs(o1[i]-o2[i]);
		dmax=(d>dmax)?d:dmax;
		v=fabs(o2[i]);
		vmax=(v>vmax)?v:vmax;
	}
	return dmax/vmax;
}

#if PERF
int
mflops (double us) {
	int i,p;
	for (i=0,p=1; i<P.d; p*=P.n[i++]);
#if DFT||INVDFT
#define C 5
#else
#define C 2.5
#endif
	return floor(C*p*log2(p)/us);
}

int
perf (void(*doit)(int)) {
	int T;
	struct timeval t1,t2;
	double us;
	int p=0;
	for (T=1; T<T_MAX; T*=2) {
		gettimeofday(&t1,NULL);
		(*doit)(T);
		gettimeofday(&t2,NULL);
		us=(t2.tv_sec-t1.tv_sec)*1000000.0+(t2.tv_usec-t1.tv_usec);
		if (us>US_MIN) {
			p=mflops(us/T);
			break;
		}
	}
	return p;
}
#endif

int
main (int argc, char **argv) {
	double r;
	if (argc!=2)
		return 1;
	if (parse(argv[1])!=0)
		return 2;
	alloc();
	init();
	setup_minfft();
	setup_fftw();
	r=cmp();
#if PERF
	int p1,p2;
	p1=perf(doit_minfft);
	p2=perf(doit_fftw);
	printf("%g %d %d\n",r,p1,p2);
#else
	printf("%g\n",r);
#endif
	cleanup_minfft();
	cleanup_fftw();
	dealloc();
	return r>THR;
}
