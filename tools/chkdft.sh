#!/bin/sh

# Run tests

# init
if [ $# -ne 1 -a $# -ne 2 ]; then
	echo "Usage: $0 test [SINGLE|EXTENDED]"
	exit 1
fi
t=$1
if [ "$2" = SINGLE ]; then
	FFTW_SFX=f
elif [ "$2" = EXTENDED ]; then
	FFTW_SFX=l
fi

# set compiler flags
# common
cflags="-std=c99 -pedantic -Wall -Wextra -lm"
# for minfft
cflags_minfft="-DMINFFT_$2 -I.. ./minfft.o"
# for FFTW
cflags_fftw="$(echo -I ~/build/fftw-mi/include -L ~/build/fftw-mi/lib -DFFTW_SFX=$FFTW_SFX -lfftw3$FFTW_SFX)"

# compile and run tests
do_tests() {
	t=$1
	for d in $2; do
		for x in $3; do
			cc -D$t -D$d -D$x chkdft.c $cflags $cflags_minfft $4
			if [ $? -ne 0 ]; then
				exit
			fi
			echo "# $t $d $x"
			./a.out
		done
	done
	rm ./a.out
}

# build minfft
echo Building minfft...
cc $cflags -DMINFFT_$2 -c ../minfft.c
echo done.

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
	gfortran -fsyntax-only -DMINFFT_$2 ../minfft.F03
	gfortran -fsyntax-only -I ~/build/fftw-mi/include fftw.f03
	# Run tests
	dims="D1 D2 D3"
	xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
	cppflags="-DFFTW_SFX=$FFTW_SFX"
	for d in $dims; do
		for x in $xforms; do
			cpp -P -D$d -D$x $cppflags chkdft.F95 -o chkdft.f95 \
			&& \
			gfortran chkdft.f95 $cflags_minfft $cflags_fftw
			if [ $? -ne 0 ]; then
				exit
			fi
			echo "# $d $x"
			./a.out
		done
	done
	rm ./a.out chkdft.f95 minfft.mod fftw.mod
fi
