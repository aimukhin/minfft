#!/bin/sh

# Run tests

# init
if [ $# -ne 1 -a $# -ne 2 ]; then
	echo "Usage: $0 test [fftw-suffix]"
	exit 1
fi
t=$1
FFTW_SFX=$2

# set compiler flags
# common
cflags="-std=c99 -pedantic -Wall -Wextra -lm -lminfft -I.. -L.."
# for FFTW
cflags_fftw="-DFFTW_SFX=$FFTW_SFX -lfftw3$FFTW_SFX"

# compile and run tests
do_tests() {
	t=$1
	for d in $2; do
		for x in $3; do
			cc -D$t -D$d -D$x chkdft.c $cflags $4
			if [ $? -ne 0 ]; then
				exit
			fi
			echo "# $t $d $x"
			./a.out
		done
	done
	rm ./a.out
}

# tests of the library itself
## inversion
if [ $t = INV ]; then
	dims="D1 D2 D3"
	xforms="DFT REALDFT DCT2 DCT3 DCT4 DST4"
	do_tests INV "$dims" "$xforms"
fi
## performance
if [ $t = PERF ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests PERF "$dims" "$xforms"
fi

# comparisons with FFTW
## accuracy
if [ $t = CMP_FFTW ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests CMP_FFTW "$dims" "$xforms" "$cflags_fftw"
fi
## performance
if [ $t = PERF_FFTW ]; then
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	do_tests PERF_FFTW "$dims" "$xforms" "$cflags_fftw"
fi

# test of Fortran interface (by comparison with FFTW)
if [ $t = FORTRAN ]; then
	# Compile interface modules
	gfortran -fsyntax-only ../minfft.f03
	gfortran -fsyntax-only fftw3.f03
	# Run tests
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	cppflags="-DFFTW_SFX=$FFTW_SFX"
	fflags="-lfftw3$FFTW_SFX"
	for d in $dims; do
		for x in $xforms; do
			cpp -P -D$d -D$x $cppflags chkdft.F95 -o chkdft.f95 \
			&& \
			gfortran chkdft.f95 -L.. -lminfft $fflags
			if [ $? -ne 0 ]; then
				exit
			fi
			echo "# $d $x"
			./a.out
		done
	done
	rm ./a.out chkdft.f95 minfft.mod fftw3.mod
fi
