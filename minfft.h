// A minimalistic FFT library
// Copyright (c) 2018 Alexander Mukhin
// MIT License

#include <complex.h>

typedef double minfft_real;
typedef double complex minfft_cmpl;

struct minfft_aux {
	int N; // number of elements to transform
	void *t; // temporary buffer
	void *e; // exponent vector
	struct minfft_aux *sub1; // subtransform structure
	struct minfft_aux *sub2; // subtransform structure
};
typedef struct minfft_aux minfft_aux;

void minfft_dft(minfft_cmpl *x, minfft_cmpl *y, const minfft_aux *a);
void minfft_invdft(minfft_cmpl *x, minfft_cmpl *y, const minfft_aux *a);
void minfft_dct2(minfft_real *x, minfft_real *y, const minfft_aux *a);
void minfft_dst2(minfft_real *x, minfft_real *y, const minfft_aux *a);
void minfft_dct3(minfft_real *x, minfft_real *y, const minfft_aux *a);
void minfft_dst3(minfft_real *x, minfft_real *y, const minfft_aux *a);
void minfft_dct4(minfft_real *x, minfft_real *y, const minfft_aux *a);
void minfft_dst4(minfft_real *x, minfft_real *y, const minfft_aux *a);

minfft_aux* minfft_mkaux_dft_1d(int N);
minfft_aux* minfft_mkaux_dft_2d(int N1, int N2);
minfft_aux* minfft_mkaux_dft_3d(int N1, int N2, int N3);
minfft_aux* minfft_mkaux_dft(int d, int *Ns);
minfft_aux* minfft_mkaux_t2t3_1d(int N);
minfft_aux* minfft_mkaux_t2t3_2d(int N1, int N2);
minfft_aux* minfft_mkaux_t2t3_3d(int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t2t3(int d, int *Ns);
minfft_aux* minfft_mkaux_t4_1d(int N);
minfft_aux* minfft_mkaux_t4_2d(int N1, int N2);
minfft_aux* minfft_mkaux_t4_3d(int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t4(int d, int *Ns);

void minfft_free_aux(minfft_aux *a);
