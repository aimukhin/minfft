// A minimalistic FFT library
// Copyright (c) 2018 Alexander Mukhin
// MIT License

#include "minfft.h"
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

extern const minfft_real pi;

// make aux data for one-dimensional forward or inverse complex DFT

//
//
//

minfft_aux*
minfft_mkaux_dft_1d (int N) {
	minfft_aux *a;
	int n;
	minfft_real *e;
	minfft_cmpl e1,e1n,e3,e3n;
	if (N<=0 || N&(N-1))
		// error if N is negative or not a power of two
		return NULL;
	a = malloc(sizeof(minfft_aux));
	if (a==NULL)
		goto err;
	a->N = N;
	if (N>=16) {
		a->t = malloc(N*sizeof(minfft_cmpl));
		if (a->t==NULL)
			goto err;
		a->e = malloc(N*sizeof(minfft_cmpl));
		if (a->e==NULL)
			goto err;
		e = a->e;
		while (N>=16) {
			for (n=0; n<N/4; n+=2) {
				e1 = exp(-2*pi*I*n/N);
				e1n = exp(-2*pi*I*(n+1)/N);
				e3 = exp(-2*pi*I*3*n/N);
				e3n = exp(-2*pi*I*3*(n+1)/N);
				*e++ = creal(e1);
				*e++ = creal(e1n);
				*e++ = creal(e3);
				*e++ = creal(e3n);
				*e++ = cimag(e1);
				*e++ = cimag(e1n);
				*e++ = cimag(e3);
				*e++ = cimag(e3n);
			}
			N /= 2;
		}
	} else {
		a->t = NULL;
		a->e = NULL;
	}
	a->sub1 = a->sub2 = NULL;
	return a;
err:	// memory allocation error
	minfft_free_aux(a);
	return NULL;
}
