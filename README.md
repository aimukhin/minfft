# fourier
This is a library of routines for computing most widely used
1-dimensional discrete Fourier transforms:

* Complex DFT and its inverse,
* DFT of real values and its inverse,
* Real symmetric transforms (DCT and DST) of types 2, 3, and 4.

It combines high performance with simplicity of implementation.

## Contents
- [Interface](#interface)
- [Transforms](#transforms)
  - [Complex DFT](#complex-dft)
  - [Inverse complex DFT](#inverse-complex-dft)
  - [Real DFT](#real-dft)
  - [Inverse real DFT](#inverse-real-dft)
  - [DCT-2](#dct-2)
  - [DST-2](#dst-2)
  - [DCT-3](#dct-3)
  - [DST-3](#dst-3)
  - [DCT-4](#dct-4)
  - [DST-4](#dst-4)
- [Implementation details](#implementation-details)
- [Performance](#performance)
- [Precision](#precision)
- [Compliance](#compliance)
- [Manually optimized versions](#manually-optimized-versions)
- [License](#license)

## Interface
All transform functions take four arguments:

* transform length `N` (a power of 2),
* pointer to the input vector `x` (its values will be destroyed),
* pointer to the output vector `y`,
* pointer to the constant *exponent vector* `e`.

The exponent vector contains precomputed constants for the given
transform type and length. The library provides functions for allocating
and filling the exponent vectors for each transform type.

All transforms are implemented in two ways - with single and double
precision. Files `dft.h` and `dft.c` contain double precision versions,
whereas `dftf.h` and `dftf.c` provide their single-precision
counterparts.

## Transforms
Here goes a description of every transform function (along with its
exponent vector generator) with a short explanation of what it
computes.

### Complex DFT
```C
void dft (int N, complex *x, complex *y, const complex *e);
complex* mkexp_dft (int N);

void dftf (int N, complex float *x, complex float *y, const complex float *e);
complex float* mkexp_dftf (int N);
```
The function computes the complex DFT, defined as usual:

![](docs/dft-def.svg)

### Inverse complex DFT
```C
void idft (int N, complex *x, complex *y, const complex *e);
complex* mkexp_idft (int N);

void idftf (int N, complex float *x, complex float *y, const complex float *e);
complex float* mkexp_idftf (int N);
```
It computes an unnormalized inverse of the complex DFT (that is, N times
the inverse).

### Real DFT
```C
void realdft (int N, double *x, double *y, const complex *e);
complex* mkexp_realdft (int N);

void realdftf (int N, float *x, float *y, const complex float *e);
complex float* mkexp_realdftf (int N);
```
This function computes the complex DFT of N reals. Since only half of N
complex outputs are independent, only those values are returned, packed
as follows:

![](docs/realdft-format.svg)

### Inverse real DFT
```C
void irealdft (int N, double *x, double *y, const complex *e);
complex* mkexp_irealdft (int N);

void irealdftf (int N, float *x, float *y, const complex float *e);
complex float* mkexp_irealdftf (int N);
```
It's an unnormalized inverse of the real DFT (N times the inverse). It
expects its input in the same format as `realdft` produces output.

### Real symmetric transforms
These transforms are actually real DFTs of a source vector extended by
symmetry to 4 or 8 times its length. Instead of writing down cumbersome
and unenlightening formulas, we'll show an example of what these
transforms do with a four-element input vector. Everywhere below, `abcd`
is the input vector `x`, and `ABCD` is the output vector `y`. Entries
marked with dots are set to zero.

#### DCT-2
```C
void dct2 (int N, double *x, double *y, const complex *e);
complex* mkexp_t2 (int N);

void dct2f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t2f (int N);
```
![](docs/dct2-def.svg)

#### DST-2
```C
void dst2 (int N, double *x, double *y, const complex *e);
complex* mkexp_t2 (int N);

void dst2f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t2f (int N);
```
![](docs/dst2-def.svg)

*NB: see [compatibility notes](#notes-on-compatibility).*

#### DCT-3
```C
void dct3 (int N, double *x, double *y, const complex *e);
complex* mkexp_t3 (int N);

void dct3f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t3f (int N);
```
![](docs/dct3-def.svg)

*NB: see [compatibility notes](#notes-on-compatibility).*

#### DST-3
```C
void dst3 (int N, double *x, double *y, const complex *e);
complex* mkexp_t3 (int N);

void dst3f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t3f (int N);
```
![](docs/dst3-def.svg)

*NB: see [compatibility notes](#notes-on-compatibility).*

#### DCT-4
```C
void dct4 (int N, double *x, double *y, const complex *e);
complex* mkexp_t4 (int N);

void dct4f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t4f (int N);
```
![](docs/dct4-def.svg)

*NB: see [compatibility notes](#notes-on-compatibility).*

#### DST-4
```C
void dst4 (int N, double *x, double *y, const complex *e);
complex* mkexp_t4 (int N);

void dst4f (int N, float *x, float *y, const complex float *e);
complex float* mkexp_t4f (int N);
```
![](docs/dst4-def.svg)

*NB: see [compatibility notes](#notes-on-compatibility).*

#### Notes on compatibility
Our definitions of real symmetric transforms differ somewhat from the
commonly used ones. Here is the summary of differences between our
definitions and those used in FFTW:

Transform | Difference
----------|-----------
`DCT-2`   | No difference
`DST-2`   | Reverse order of outputs, our transform differs by constant `-1`
`DCT-3`   | Our transform differs by constant `2`
`DST-3`   | Reverse order of inputs, our transform differs by constant `-2`
`DCT-4`   | Our transform differs by constant `2`
`DST-4`   | Our transform differs by constant `-2`

## Implementation details
The complex DFT is computed by a split-radix (2/4)
decimation-in-frequency explicitly recursive fast Fourier transform.
This method achieves a remarkable balance between performance and
simplicity, and it behaves particularly cache-friendly, since it refers
mostly to adjacent memory locations.

All the real transforms are reduced eventually to a half-length complex
transform. Further details are given [here](docs/math-details.md).

## Performance
These graphs show the execution times of our routines compared with
those of FFTW. Here our library is compiled with GCC, which was unable
to completely vectorize the code, whereas FFTW codelets are vectorized.
The details of the test environment are given
[below](#test-environment).

Timings are measured in microseconds per single call, and are plotted
with three-sigma error intervals.

### Complex DFT
![](docs/perf-dft.svg)

### Real DFT
![](docs/perf-realdft.svg)

For the complex DFTs, our library performs slower than FFTW by 1.5-2
times (and even more for the single-precision routines). This is not
surprising, considering the sophistication of FFTW and lack of
vectorization ability for this particular compiler.

Nevertheless, for the real transforms things begin to look better, and
our performance approaches to, and sometimes exceeds that of FFTW. One
can notice two areas where our library performs better - for the very
small and very big transforms.

### DCT-2
![](docs/perf-dct2.svg)

### DCT-4
![](docs/perf-dct4.svg)

For the real symmetric transforms, our library performs at least as well
as FFTW for most transform sizes, and sometimes much better.

How can that be, provided that our core complex DFT performs worse than
FFTW? Probably this is because we use more efficient ways of reducing
these transforms to the core complex one.

## Precision
This graph shows how much the results of real transforms differ
from the results obtained by direct computation of the underlying
complex DFT. Absent data mean the difference is exactly zero.

![](docs/prec.svg)

And this graph shows the distribution of an absolute error within the
transform results. Here, the DCT-4 of length N=1024 is compared with its
direct computation by the complex DFT of length 8N:

![](docs/prec-dct4-1024.svg)

The other transforms exhibit the similar uniform error distribution.

## Test environment
The above tests and comparisons are made with the current library source
compiled with GCC version 7.3.0 on x86_64 with the only optimization
option `-Ofast`. Unfortunately, the compiler complained that it was
unable to vectorize some parts of the code, and sometimes decided that
it would not be profitable.

The version of FFTW used is 3.3.7-1 packed for Arch Linux x86_64. FFTW
plans are created with options FFTW_ESTIMATE and FFTW_DESTROY_INPUT.

The performance measurements are made on an isolated core of an Intel®
Celeron® N3050 CPU running at 2160 MHz.

The source data files for the graphs, along with the programs `chkdft.c`
and `chkdftf.c` used to conduct the above (and many other) tests, are
available in the `tests` subdirectory.

The `chkdft.c` program contains a lot of examples of library functions
usage, and therefore can serve as a reference.

## Compliance
The source code complies with the C99 standard.

## Manually optimized versions
We also provide manually written assembly-language implementations of
the forward and inverse complex DFTs. Since they are the core transforms
to which the other transforms are ultimately reduced, it's worth to
invest some effort in computing them as fast as possible.

The performance of the hand-written code is much better than of the code
emitted by GCC.

Machine-dependent versions are kept in separate branches. At present,
we provide the following version:

* For [x86-64 with SSE3](../../tree/x86-64-sse3-sysv), using SystemV ABI calling conventions.

## License
The library is in the public domain.
