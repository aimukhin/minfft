#!/bin/sh

# Compiler flags
cflags="-std=c99 -pedantic -Wall -Wextra -lm -lminfft -I.. -L.."
cflags_fftw="-lfftw3" # double
#cflags_fftw="-DFFTW_PFX=fftwf -lfftw3f" # single
#cflags_fftw="-DFFTW_PFX=fftwl -lfftw3l" # long double
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
	echo "Usage: $0 test|clean"
	exit 1
fi
t=$1

if [ $t = clean ]; then
	rm -f a.out *.mod chkdft.f95
fi

# Tests of the library itself
## Inversion
if [ $t = INV ]; then
	dims="D1 D2 D3"
	xforms="DFT REALDFT DCT2 DCT3 DCT4 DST4"
	do_tests INV "$dims" "$xforms"
fi
## Performance
if [ $t = PERF ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests PERF "$dims" "$xforms"
fi

# Comparisons with FFTW
## Accuracy
if [ $t = CMP_FFTW ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests CMP_FFTW "$dims" "$xforms" "$cflags_fftw"
fi
## Performance
if [ $t = PERF_FFTW ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests PERF_FFTW "$dims" "$xforms" "$cflags_fftw"
fi

# Comparisons with Kiss FFT
## Accuracy
if [ $t = CMP_KISS ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests CMP_KISS "$dims" "$xforms" "$cflags_kiss"
fi
## Performance
if [ $t = PERF_KISS ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests PERF_KISS "$dims" "$xforms" "$cflags_kiss"
fi

# Comparisons with Ne10
## Accuracy
if [ $t = CMP_NE10 ]; then
	dims="D1"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests CMP_NE10 "$dims" "$xforms" "$cflags_ne10"
fi
## Performance
if [ $t = PERF_NE10 ]; then
	dims="D1"
	xforms="DFT INVDFT REALDFT INVREALDFT"
	do_tests PERF_NE10 "$dims" "$xforms" "$cflags_ne10"
fi

# Test of Fortran interface (by comparison with FFTW)
if [ $t = FORTRAN ]; then
	# Compile interface modules
	gfortran -fsyntax-only ../minfft.f03
	gfortran -fsyntax-only fftw3.f03
	# Run tests
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
# double precision
	fflags="-lfftw3"
# single precision
#	cppflags="-DFFTW_PFX=fftwf"
#	fflags="-lfftw3f"
# extended precision
#	cppflags="-DFFTW_PFX=fftwl"
#	fflags="-lfftw3l"
	for d in $dims; do
		for x in $xforms; do
			cpp -P -D$d -D$x $cppflags chkdft.F95 -o chkdft.f95
			gfortran chkdft.f95 -L.. -lminfft $fflags
			echo "# $d $x"
			./a.out
		done
	done
fi
