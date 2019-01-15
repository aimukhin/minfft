# minfft
A minimalistic Fast Fourier Transform library.

It achieves high performance by simple means.

## Overview

The library routines compute:

* Forward and inverse complex DFT,
* Forward and inverse DFT of real data,
* Cosine and sine transforms of the types 2, 3, 4

of any dimensionality and power-of-two lengths.

## Contents
- [Interface](#interface)
- [Data types](#data-types)
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
- [Freeing auxiliary data](#freeing-auxiliary-data)
- [Memory requirements](#memory-requirements)
- [Implementation details](#implementation-details)
- [Performance](#performance)
- [Test environment](#test-environment)
- [Conformance](#conformance)
- [SIMD branches](#simd-branches)
- [License](#license)

## Interface
All transform routines take three arguments:

* a pointer to the input data `x`,
* a pointer to the output data `y`,
* a pointer to the auxiliary data `a`.

The transform routines are capable of both in-place and out-of-place
operation. In the latter case the input data would be left intact.

Auxiliary data contain chains of precomputed constants and temporary
memory buffers, required by a transform routine to do its job. Once
prepared, the auxiliary data can be reused as many times as needed.
Also, the same auxiliary data fit for both forward and inverse
transforms of the same kind.

Here is an example to give you an idea how the library functions are
used:
```C
	double complex x[N],y[N];
	// prepare aux data
	minfft_aux *a = minfft_mkaux_dft_1d(N);
	// do transforms
	minfft_dft(x,y,a);
	minfft_invdft(y,x,a);
	// free aux data
	minfft_free_aux(a);
```

## Data types
By default, the library operates with double precision:
```C
	typedef double minfft_real;
	typedef double complex minfft_cmpl;
```
By changing these definitions in the header file, you can adapt the
library to use `float` or `long double` instead.

## Transforms
Below is a list of transform functions and their auxiliary data
makers. For convenience, we provide makers for one-, two- and
three-dimensional transforms, along with a generic any-dimensional one.
The order of dimensions follows the C way — the last index changes most
quickly. Auxiliary data makers return NULL if an error occured.

Also we give a formal definition of each transform for the
one-dimensional case.

Our definitions of transforms, and input and output data format of the
transform routines, are fully compatible with FFTW.

### Complex DFT
![](docs/dft.svg)
```C
minfft_aux* minfft_mkaux_dft_1d (int N);
minfft_aux* minfft_mkaux_dft_2d (int N1, int N2);
minfft_aux* minfft_mkaux_dft_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_dft (int d, int *Ns);
void minfft_dft (minfft_cmpl *x, minfft_cmpl *y, const minfft_aux *a);
```

### Inverse complex DFT
![](docs/invdft.svg)
```C
minfft_aux* minfft_mkaux_dft_1d (int N);
minfft_aux* minfft_mkaux_dft_2d (int N1, int N2);
minfft_aux* minfft_mkaux_dft_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_dft (int d, int *Ns);
void minfft_invdft (minfft_cmpl *x, minfft_cmpl *y, const minfft_aux *a);
```

### Real DFT
This transform returns mostly non-redundant part of the complex DFT of
real data.

For a real array of dimensions

![](docs/realdft-in.svg)

it produces a complex array of dimensions

![](docs/realdft-out.svg)

Note that output takes a little more space than input. For in-place
operation, make sure the data buffer is big enough to contain output.

```C
minfft_aux* minfft_mkaux_realdft_1d (int N);
minfft_aux* minfft_mkaux_realdft_2d (int N1, int N2);
minfft_aux* minfft_mkaux_realdft_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_realdft (int d, int *Ns);
void minfft_realdft (minfft_real *x, minfft_cmpl *z, const minfft_aux *a);
```

### Inverse real DFT
This is the inversion of the real DFT.

It takes a complex array of dimensions 

![](docs/realdft-out.svg)

and returns a real array of dimensions

![](docs/realdft-in.svg)

It is better fed with the output of the real DFT routine, especially in
multi-dimensional case, since it expects input with some residual
redundancies.

```C
minfft_aux* minfft_mkaux_realdft_1d (int N);
minfft_aux* minfft_mkaux_realdft_2d (int N1, int N2);
minfft_aux* minfft_mkaux_realdft_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_realdft (int d, int *Ns);
void minfft_invrealdft (minfft_cmpl *z, minfft_real *y, const minfft_aux *a);
```

#### DCT-2
![](docs/dct2.svg)
```C
minfft_aux* minfft_mkaux_t2t3_1d (int N);
minfft_aux* minfft_mkaux_t2t3_2d (int N1, int N2);
minfft_aux* minfft_mkaux_t2t3_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t2t3 (int d, int *Ns);
void minfft_dct2 (minfft_real *x, minfft_real *y, const minfft_aux *a);
```

#### DST-2
![](docs/dst2.svg)
```C
minfft_aux* minfft_mkaux_t2t3_1d (int N);
minfft_aux* minfft_mkaux_t2t3_2d (int N1, int N2);
minfft_aux* minfft_mkaux_t2t3_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t2t3 (int d, int *Ns);
void minfft_dst2 (minfft_real *x, minfft_real *y, const minfft_aux *a);
```

#### DCT-3
![](docs/dct3.svg)
```C
minfft_aux* minfft_mkaux_t2t3_1d (int N);
minfft_aux* minfft_mkaux_t2t3_2d (int N1, int N2);
minfft_aux* minfft_mkaux_t2t3_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t2t3 (int d, int *Ns);
void minfft_dct3 (minfft_real *x, minfft_real *y, const minfft_aux *a);
```

#### DST-3
![](docs/dst3.svg)
```C
minfft_aux* minfft_mkaux_t2t3_1d (int N);
minfft_aux* minfft_mkaux_t2t3_2d (int N1, int N2);
minfft_aux* minfft_mkaux_t2t3_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t2t3 (int d, int *Ns);
void minfft_dst3 (minfft_real *x, minfft_real *y, const minfft_aux *a);
```

#### DCT-4
![](docs/dct4.svg)
```C
minfft_aux* minfft_mkaux_t4_1d (int N);
minfft_aux* minfft_mkaux_t4_2d (int N1, int N2);
minfft_aux* minfft_mkaux_t4_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t4 (int d, int *Ns);
void minfft_dct4 (minfft_real *x, minfft_real *y, const minfft_aux *a);
```

#### DST-4
![](docs/dst4.svg)
```C
minfft_aux* minfft_mkaux_t4_1d (int N);
minfft_aux* minfft_mkaux_t4_2d (int N1, int N2);
minfft_aux* minfft_mkaux_t4_3d (int N1, int N2, int N3);
minfft_aux* minfft_mkaux_t4 (int d, int *Ns);
void minfft_dst4 (minfft_real *x, minfft_real *y, const minfft_aux *a);
```

## Freeing auxiliary data
If not needed anymore, the memory consumed by the auxiliary data can
be freed by the `minfft_free_aux()` routine:
```C
void minfft_free_aux (minfft_aux *a);
```

## Memory requirements
Our library does not try to save memory, and allocates temporary buffers
where it benefits performance.

The amounts of memory, allocated for the auxiliary data of the
one-dimensional transforms, are given below:

Transform                                | Auxiliary data size
-----------------------------------------|---------------------
Complex DFT of length `N`                | `2N` complex numbers
Real DFT of length `N`                   | `3.5N` real numbers
Type-2 or Type-3 transform of length `N` | `5.5N` real numbers
Type-4 transform of length `N`           | `6N` real numbers

Multi-dimensional transforms use a temporary buffer of the same size as
the input data. This value is the dominant term in their auxiliary data
size.

## Implementation details
The complex DFT is computed by a split-radix (2/4), decimation in
frequency, explicitly recursive fast Fourier transform. This method
achieves a remarkable balance between performance and simplicity, and it
behaves particularly cache-friendly, since it refers mostly to adjacent
memory locations.

The real transforms are reduced eventually to a half-length complex
transform.

For each transform, we first implement its one-dimensional,
out-of-place, input-preserving, sequential input, strided output
routine. This allows us to compute a multi-dimensional transform by
repeated application of its one-dimensional routine along each
dimension.

## Performance
Below is a plot of the execution times of our one-dimensional complex
DFT routine, compared with FFTW and Kiss FFT. The results are measured
in microseconds per call.

![](docs/perf-dft-1d.svg)

## Test environment
The libraries being compared are built with the GNU C compiler version
8.2.1 for the x86_64 platform. The only optimization option set is
`-Ofast`.

The version of FFTW used is 3.3.8. To make a fair comparison, we disable
all its SIMD optimizations, and therefore compare performance of the
machine-independent code. FFTW plans are created with the option
FFTW_ESTIMATE.

The version of Kiss FFT used is 1.3.0.

The performance measurements are made on an isolated core of an Intel®
Celeron® N3050 CPU running at 2160 MHz.

The results of many performance and accuracy tests, along with the
programs used to conduct them, are available in the `tests` subdirectory.

## Conformance
The source code conforms to the C99 standard.

## SIMD branches
We also provide machine-dependent branches, optimized manually for SIMD
architectures. These are:

- [neon-single](../../tree/neon-single)
- [sse3-double](../../tree/sse3-double)
- [sse3-single](../../tree/sse3-single)

## License
MIT.
