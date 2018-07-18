// DFT library (double-precision version).
// Written by Alexander Mukhin.
// Public domain.

#include <complex.h>

struct aux {
	int N; // number of elements to transform
	void *t; // temporary buffer
	void *e; // exponent vector
	struct aux *sub1; // next structure
	struct aux *sub2; // next structure
};
void dft(double complex *x, double complex *y, const struct aux *a);
void idft(double complex *x, double complex *y, const struct aux *a);
void realdft_1d(double *x, double *y, const struct aux *a);
void irealdft_1d(double *x, double *y, const struct aux *a);
void dct2(double *x, double *y, const struct aux *a);
void dst2(double *x, double *y, const struct aux *a);
void dct3(double *x, double *y, const struct aux *a);
void dst3(double *x, double *y, const struct aux *a);
void dct4(double *x, double *y, const struct aux *a);
void dst4(double *x, double *y, const struct aux *a);
struct aux *mkaux_dft(int d, int *n);
struct aux *mkaux_idft(int d, int *n);
struct aux *mkaux_realdft_1d(int N);
struct aux *mkaux_irealdft_1d(int N);
struct aux *mkaux_t2(int d, int *n);
struct aux *mkaux_t3(int d, int *n);
struct aux *mkaux_t4(int d, int *n);
void free_aux(struct aux *a);
