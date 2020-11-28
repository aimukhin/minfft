#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>

#define D_MAX 3
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
	int d;
	int n[D_MAX];
	void *in;
	void *out;
	void *x;
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
	P.in=malloc(p*sizeof(minfft_cmpl));
	P.out=malloc(p*sizeof(minfft_cmpl));
#elif REALDFT
	for (p=1,i=0; i<P.d-1; p*=P.n[i++]);
	P.in=malloc(p*P.n[P.d-1]*sizeof(minfft_real));
	P.out=malloc(p*(P.n[P.d-1]/2+1)*sizeof(minfft_cmpl));
#elif INVREALDFT
	for (p=1,i=0; i<P.d-1; p*=P.n[i++]);
	P.in=malloc(p*(P.n[P.d-1]/2+1)*sizeof(minfft_cmpl));
	P.out=malloc(p*P.n[P.d-1]*sizeof(minfft_real));
#elif DCT2||DST2||DCT3||DST3||DCT4||DST4
	for (p=1,i=0; i<P.d; p*=P.n[i++]);
	P.in=malloc(p*sizeof(minfft_real));
	P.out=malloc(p*sizeof(minfft_real));
#endif
}

void
setup_minfft (void) {
#if DFT||INVDFT
	P.x=minfft_mkaux_dft(P.d,P.n);
#elif REALDFT||INVREALDFT
	P.x=minfft_mkaux_realdft(P.d,P.n);
#elif DCT2||DST2||DCT3||DST3
	P.x=minfft_mkaux_t2t3(P.d,P.n);
#elif DCT4||DST4
	P.x=minfft_mkaux_t4(P.d,P.n);
#endif
}

void
doit_minfft (int T) {
	int t;
	for (t=0; t<T; ++t)
#if DFT
		minfft_dft(P.in,P.out,P.x);
#elif INVDFT
		minfft_invdft(P.in,P.out,P.x);
#elif REALDFT
		minfft_realdft(P.in,P.out,P.x);
#elif INVREALDFT
		minfft_invrealdft(P.in,P.out,P.x);
#elif DCT2
		minfft_dct2(P.in,P.out,P.x);
#elif DST2
		minfft_dst2(P.in,P.out,P.x);
#elif DCT3
		minfft_dct3(P.in,P.out,P.x);
#elif DST3
		minfft_dst3(P.in,P.out,P.x);
#elif DCT4
		minfft_dct4(P.in,P.out,P.x);
#elif DST4
		minfft_dst4(P.in,P.out,P.x);
#endif
}

void
cleanup_minfft (void) {
	minfft_free_aux(P.x);
}

void
setup_fftw (void) {
#if DFT
	P.x=FFTW(plan_dft)(P.d,P.n,P.in,P.out,FFTW_FORWARD,FFTW_ESTIMATE);
#elif INVDFT
	P.x=FFTW(plan_dft)(P.d,P.n,P.in,P.out,FFTW_BACKWARD,FFTW_ESTIMATE);
#elif REALDFT
	P.x=FFTW(plan_dft_r2c)(P.d,P.n,P.in,P.out,FFTW_ESTIMATE);
#elif INVREALDFT
	P.x=FFTW(plan_dft_c2r)(P.d,P.n,P.in,P.out,FFTW_ESTIMATE);
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
	P.x=FFTW(plan_r2r)(P.d,P.n,P.in,P.out,kinds,FFTW_ESTIMATE);
#endif
}

void
doit_fftw (int T) {
	int t;
	for (t=0; t<T; ++t)
		FFTW(execute)(P.x);
}

void
cleanup_fftw (void) {
	FFTW(destroy_plan)(P.x);
}

void
dealloc (void) {
	free(P.in);
	free(P.out);
}

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
perf (void(*setup)(void), void(*doit)(int), void(*cleanup)(void)) {
	int T;
	struct timeval t1,t2;
	double us;
	int p=-1;
	(*setup)();
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
	(*cleanup)();
	return p;
}

int
main (int argc, char **argv) {
	int p1,p2;
	if (argc!=2)
		return 1;
	if (parse(argv[1])!=0)
		return 2;
	alloc();
	p1=perf(setup_minfft,doit_minfft,cleanup_minfft);
	p2=perf(setup_fftw,doit_fftw,cleanup_fftw);
	dealloc();
	if (p1<0||p2<0)
		return 3;
	printf("%d %d\n",p1,p2);
	return 0;
}
