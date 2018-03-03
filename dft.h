// DFT library
// Written by Alexander Mukhin
// Public domain

#include <complex.h>

void dft (int N, complex *x, complex *y, const complex *e);
complex* mkexp_dft (int N);

void idft (int N, complex *x, complex *y, const complex *e);
complex* mkexp_idft (int N);

void realdft (int N, double *x, double *y, const complex *e);
complex* mkexp_realdft (int N);

void irealdft (int N, double *x, double *y, const complex *e);
complex* mkexp_irealdft (int N);

void dct2 (int N, double *x, double *y, const complex *e);
void dst2 (int N, double *x, double *y, const complex *e);
complex* mkexp_t2 (int N);

void dct3 (int N, double *x, double *y, const complex *e);
void dst3 (int N, double *x, double *y, const complex *e);
complex* mkexp_t3 (int N);

void dct4 (int N, double *x, double *y, const complex *e);
void dst4 (int N, double *x, double *y, const complex *e);
complex* mkexp_t4 (int N);
