
#include "minfft.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))

/* random number generator, 0 <= RND < 1 */
#define RND(p) ((*(p) = (*(p) * 7141 + 54773) % 259200) * (1.0 / 259200.0))

#if MINFFT_SINGLE
#define MINFFT_PREC_TYPE  "single"
#define ERR_LIMIT 2.0e-6
#define MIN_FABS(x)  fabsf(x)
#define MIN_COS(x)   cosf(x)
#define MIN_SIN(x)   sinf(x)

#elif MINFFT_EXTENDED
#define MINFFT_PREC_TYPE  "long double"
#define ERR_LIMIT 3.0e-15
#define MIN_FABS(x)  fabsl(x)
#define MIN_COS(x)   cosl(x)
#define MIN_SIN(x)   sinl(x)

#else
#define MINFFT_PREC_TYPE  "double"
#define ERR_LIMIT 3.0e-15
#define MIN_FABS(x)  fabs(x)
#define MIN_COS(x)   cos(x)
#define MIN_SIN(x)   sin(x)

#endif


// constants
static const minfft_real pi=3.141592653589793238462643383279502884L;

void putdata(int N, minfft_real *a);
void printdata(int N, minfft_real scale, const minfft_real *a, int print_cols);
minfft_real errorcheck(int N, minfft_real scale, const minfft_real *a, int *maxindex);
minfft_real abs_max(int N, const minfft_real *a);
int max_idx(int N, const minfft_real *a);


static inline minfft_real * c2r(minfft_cmpl * c)
{
    void *v = c;
    minfft_real * r = (minfft_real *)c;
    return r;
}


int main(int argc, char *argv[])
{
    int k, j, n, ret, retCode;
    int NMAX, max_err_index;
    minfft_real *xr, *xi;
    minfft_cmpl *Xr, *zr, *Xz, *zi;
    minfft_aux *prep_min_ftz, *prep_min_ftr;
    minfft_real err, mx;

    printf("minfft %s-precision FFT-library and sizeof(minfft_real) = %u\n", MINFFT_PREC_TYPE, (unsigned)sizeof(minfft_real));

    if (1 < argc)
    {
        n = atoi(argv[1]);
        printf("running %s with N = %d\n", argv[0], n);
    }
    else
    {
        fprintf(stderr, "usage: %s <N>\n", argv[0]);
        fprintf(stderr, "  N = data length: must be 2^m\n");
        return 1;
    }

    NMAX = n;

    /* xr: real input; Xr: spectrum of xr; xi: inverse from spectrum */
    xr = (minfft_real*)malloc( NMAX * sizeof(minfft_real));
    Xr = (minfft_cmpl*)malloc( (NMAX / 2 + 1) * sizeof(minfft_cmpl));
    xi = (minfft_real*)malloc( NMAX * sizeof(minfft_real));
    /* zr: complex input; Xr: spectrum of zr; zi: inverse from spectrum */
    zr = (minfft_cmpl*)malloc( NMAX * sizeof(minfft_cmpl));
    Xz = (minfft_cmpl*)malloc( NMAX * sizeof(minfft_cmpl));
    zi = (minfft_cmpl*)malloc( NMAX * sizeof(minfft_cmpl));

    retCode = 0;

    /* check of CDFT */
    putdata(2 * n, c2r(zr));
    mx = abs_max(2 * n, c2r(zr));

    prep_min_ftz = minfft_mkaux_dft_1d(n);
    minfft_dft(zr, Xz, prep_min_ftz);
    minfft_invdft(Xz, zi, prep_min_ftz);

    /* check: inverse of forward transformed  equals input? */
    err = errorcheck(2 * n, 1.0 / n, c2r(zi), &max_err_index);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("complex dft max err= %g%s at index %d for max(abs(inp)) = %g\n",
        (double)err, ret ? "":" (within error tolerance)", max_err_index, (double)mx);
    if (ret) {
        printf("\ninput:\n");
        printdata(2 * n, 1.0, c2r(zr), 4);
        printf("\ninverse of spectrum:\n");
        printdata(2 * n, 1.0 / n, c2r(zi), 4);
    }

    /* check: forward of carrier produces single peak in spectrum at right position? */
    for (k = 0; k < n; k += MAX(1, n/8))
    {
        minfft_real *t = c2r(zr);
        minfft_real *z = c2r(Xz);
        const minfft_real dphi = 2 * pi * k / n;
        const minfft_real scale = (minfft_real)1.0 / n;
        minfft_real maxval, max2;
        const int specsize = n;
        int mi, ni;
        /* prepare input: single carrier */
        for (j = 0; j < n; ++j) {
            t[2*j] = MIN_COS(j*dphi);
            t[2*j+1] = MIN_SIN(j*dphi);
        }
        /* forward fft */
        minfft_dft(zr, Xz, prep_min_ftz);
        /* calculate power; reuse xi[] */
        for (j = 0; j < specsize; ++j)
            xi[j] = (scale * z[2*j]) * (scale * z[2*j]) + (scale * z[2*j+1]) * (scale * z[2*j+1]);
        /* find / check maxima: position and ratio */
        mi = max_idx(specsize, xi);
        maxval = xi[mi];
        xi[mi] = 0;
        ni = max_idx(specsize, xi);
        max2 = xi[ni];

        ret = (mi == k && max2 * 1000 < maxval) ? 0 : 1;
        retCode += ret;
        printf("complex carrier %d: max %g at %d;  max2 %g at %d --> err %d\n",
            k, (double)maxval, mi, (double)max2, ni, ret);
    }
    minfft_free_aux(prep_min_ftz);

    /* check of RDFT */
    putdata(n, xr);
    mx = abs_max(n, xr);

    prep_min_ftr = minfft_mkaux_realdft_1d(n);
    minfft_realdft(xr, Xr, prep_min_ftr);
    minfft_invrealdft(Xr, xi, prep_min_ftr);

    /* check: inverse of forward transformed  equals input? */
    err = errorcheck(n, 1.0 / n, xi, &max_err_index);
    ret = (err > ERR_LIMIT) ? 1 : 0;
    retCode += ret;
    printf("real dft max err= %g%s at index %d for max(abs(inp)) = %g\n",
        (double)err, ret ? "":" (within error tolerance)", max_err_index, (double)mx);
    if (ret) {
        printf("\ninput:\n");
        printdata(n, 1.0, xr, 4);
        printf("\ninverse of spectrum:\n");
        printdata(n, 1.0 / n, xi, 4);
    }

    /* check: forward of carrier produces single peak in spectrum at right position? */
    for (k = 0; k <= n/2; k += MAX(1, n/8))
    {
        minfft_real *z = c2r(Xr);
        const minfft_real dphi = 2 * pi * k / n;
        const minfft_real scale = (minfft_real)1.0 / n;
        minfft_real maxval, max2;
        const int specsize = n /2 +1;
        int mi, ni;
        /* prepare input: single carrier */
        for (j = 0; j < n; ++j) {
            xr[j] = MIN_COS(j*dphi);
        }
        /* forward fft */
        minfft_realdft(xr, Xr, prep_min_ftr);

        /* calculate power; reuse xi[] */
        for (j = 0; j < specsize; ++j)
            xi[j] = (scale * z[2*j]) * (scale * z[2*j]) + (scale * z[2*j+1]) * (scale * z[2*j+1]);
        /* find / check maxima: position and ratio */
        mi = max_idx(specsize, xi);
        maxval = xi[mi];
        xi[mi] = 0;
        ni = max_idx(specsize, xi);
        max2 = xi[ni];

        ret = (mi == k && max2 * 1000 < maxval) ? 0 : 1;
        retCode += ret;
        printf("real carrier %d: max %g at %d;  max2 %g at %d --> err %d\n",
            k, (double)maxval, mi, (double)max2, ni, ret);
    }
    minfft_free_aux(prep_min_ftr);

    free(xr);
    free(Xr);
    free(xi);
    free(zr);
    free(Xz);
    free(zi);

    return retCode;
}


void putdata(int N, minfft_real *a)
{
    int j, seed = 0;

    for (j = 0; j < N; j++) {
        a[j] = RND(&seed);
    }
}

void printdata(int N, minfft_real scale, const minfft_real *a, int print_cols)
{
    int j;
    if (!print_cols)
        return;

    for (j = 0; j < N; j++) {
        if (!(j % print_cols))
            printf("\n%2d: ", j);
        printf("%g%s", (double)(scale * a[j]), (j % print_cols) ? "\t":" ");
    }
    printf("\n");
}

minfft_real errorcheck(int N, minfft_real scale, const minfft_real *a, int *maxindex)
{
    minfft_real err = 0, e;
    int j, seed = 0;
    *maxindex = 0;

    for (j = 0; j < N; j++) {
        e = RND(&seed) - a[j] * scale;
        if (MIN_FABS(e) > err)
            *maxindex = j;
        err = MAX(err, MIN_FABS(e));
    }
    return err;
}

minfft_real abs_max(int N, const minfft_real *a)
{
    minfft_real m = 0;
    int j = 0;

    for (j = 0; j < N; j++) {
        m = MAX(m, MIN_FABS(a[j]));
    }
    return m;
}

int max_idx(int N, const minfft_real *a)
{
    minfft_real m = 0;
    int j = 0, mi = 0;

    for (j = 0; j < N; j++) {
        if (a[j] > m) {
            m = a[j];
            mi = j;
        }
    }
    return mi;
}

