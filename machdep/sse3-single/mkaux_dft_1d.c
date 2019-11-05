// A minimalistic FFT library
// Copyright (c) 2018 Alexander Mukhin
// MIT License

#include "minfft.h"
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

extern const minfft_real pi;

// make aux data for one-dimensional forward or inverse complex DFT

// a modified routine with two iterations unrolled and values rearranged
// to facilitate parallel processing in assembly-language code

minfft_aux*
minfft_mkaux_dft_1d (int N) {
	minfft_aux *a;
	int n;
	minfft_cmpl *e;
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
				*e++ = exp(-2*pi*I*n/N);
				*e++ = exp(-2*pi*I*(n+1)/N);
				*e++ = exp(-2*pi*I*3*n/N);
				*e++ = exp(-2*pi*I*3*(n+1)/N);
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
