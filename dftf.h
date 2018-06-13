// DFT library (single-precision version).
// Written by Alexander Mukhin.
// Public domain.

#include <complex.h>

void dftf (int N, float complex *x, float complex *y, const float complex *e);
float complex* mkexp_dftf (int N);

void idftf (int N, float complex *x, float complex *y, const float complex *e);
float complex* mkexp_idftf (int N);

void realdftf (int N, float *x, float *y, const float complex *e);
float complex* mkexp_realdftf (int N);

void irealdftf (int N, float *x, float *y, const float complex *e);
float complex* mkexp_irealdftf (int N);

void dct2f (int N, float *x, float *y, const float complex *e);
void dst2f (int N, float *x, float *y, const float complex *e);
float complex* mkexp_t2f (int N);

void dct3f (int N, float *x, float *y, const float complex *e);
void dst3f (int N, float *x, float *y, const float complex *e);
float complex* mkexp_t3f (int N);

void dct4f (int N, float *x, float *y, const float complex *e);
void dst4f (int N, float *x, float *y, const float complex *e);
float complex* mkexp_t4f (int N);
