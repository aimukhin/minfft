// DFT library (single-precision version).
// Written by Alexander Mukhin.
// Public domain.

#include <complex.h>

void dftf (int N, complex float *x, complex float *y, const complex float *e);
complex float* mkexp_dftf (int N);

void idftf (int N, complex float *x, complex float *y, const complex float *e);
complex float* mkexp_idftf (int N);

void realdftf (int N, float *x, float *y, const complex float *e);
complex float* mkexp_realdftf (int N);

void irealdftf (int N, float *x, float *y, const complex float *e);
complex float* mkexp_irealdftf (int N);

void dct2f (int N, float *x, float *y, const complex float *e);
void dst2f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t2f (int N);

void dct3f (int N, float *x, float *y, const complex float *e);
void dst3f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t3f (int N);

void dct4f (int N, float *x, float *y, const complex float *e);
void dst4f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t4f (int N);
