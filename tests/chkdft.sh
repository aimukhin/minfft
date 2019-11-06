#!/bin/sh

# Compiler flags
cflags="-std=c99 -pedantic -Wall -Wextra -lm -lminfft -I.. -L.."
cflags_fftw="-lfftw3"
#cflags_fftw="-DFFTW_PFX=fftwf -lfftw3f"
cflags_kiss=$(echo -I ~/build/kiss_fft/ ~/build/kiss_fft/kiss_fft*.o)
cflags_ne10=$(echo -std=gnu17 -I ~/build/Ne10/inc/ -L ~/build/Ne10/modules/ -lNE10)

# Compile and run tests
do_tests() {
	t=$1
	for d in $2; do
		for x in $3; do
			cc -D$t -D$d -D$x chkdft.c $cflags $4
			echo "# $t $d $x"
			./a.out
		done
	done
}

if [ $# -ne 1 ]; then
	echo "Usage: $0 <test|ALL>"
	exit 1
fi
t=$1

# Tests of the library itself
## Inversion
if [ $t = INV -o $t = ALL ]; then
	dims="D1 D2 D3"
	xforms="DFT REALDFT DCT2 DCT3 DCT4 DST4"
	do_tests INV "$dims" "$xforms"
fi
## Performance
if [ $t = PERF -o $t = ALL ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests PERF "$dims" "$xforms"
fi

# Comparisons with FFTW
## Accuracy
if [ $t = CMP_FFTW -o $t = ALL ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests CMP_FFTW "$dims" "$xforms" "$cflags_fftw"
fi
## Performance
if [ $t = PERF_FFTW -o $t = ALL ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests PERF_FFTW "$dims" "$xforms" "$cflags_fftw"
fi

# Comparisons with Kiss FFT
## Accuracy
if [ $t = CMP_KISS -o $t = ALL ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests CMP_KISS "$dims" "$xforms" "$cflags_kiss"
fi
## Performance
if [ $t = PERF_KISS -o $t = ALL ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests PERF_KISS "$dims" "$xforms" "$cflags_kiss"
fi

# Comparisons with Ne10
## Accuracy
if [ $t = CMP_NE10 -o $t = ALL ]; then
	dims="D1"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests CMP_NE10 "$dims" "$xforms" "$cflags_ne10"
fi
## Performance
if [ $t = PERF_NE10 -o $t = ALL ]; then
	dims="D1"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests PERF_NE10 "$dims" "$xforms" "$cflags_ne10"
fi

# Test of Fortran interface (by comparison with FFTW)
if [ $t = FORTRAN -o $t = ALL ]; then
	# Compile interface modules
	gfortran -fsyntax-only ../minfft.f03
	gfortran -fsyntax-only fftw3.f03
	# Run tests
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	fflags="-L.. -lminfft -lfftw3"
	for d in $dims; do
		for x in $xforms; do
			gfortran -D$d -D$x chkdft.F95 $fflags
			echo "# $d $x"
			./a.out
		done
	done
fi
