// DFT library.
// Written by Alexander Mukhin.
// Public domain.

#include <complex.h>

struct dft_aux {
	int N; // number of elements to transform
	void *t; // temporary buffer
	void *e; // exponent vector
	struct dft_aux *sub1; // subtransform structure
	struct dft_aux *sub2; // subtransform structure
};
void dft(double complex *x, double complex *y, const struct dft_aux *a);
void invdft(double complex *x, double complex *y, const struct dft_aux *a);
void realdft_1d(double *x, double *y, const struct dft_aux *a);
void invrealdft_1d(double *x, double *y, const struct dft_aux *a);
void dct2(double *x, double *y, const struct dft_aux *a);
void dst2(double *x, double *y, const struct dft_aux *a);
void dct3(double *x, double *y, const struct dft_aux *a);
void dst3(double *x, double *y, const struct dft_aux *a);
void dct4(double *x, double *y, const struct dft_aux *a);
void dst4(double *x, double *y, const struct dft_aux *a);
struct dft_aux *mkaux_complex(int d, int *Ns);
struct dft_aux *mkaux_real_1d(int N);
struct dft_aux *mkaux_t2t3(int d, int *Ns);
struct dft_aux *mkaux_t4(int d, int *Ns);
void free_aux(struct dft_aux *a);
