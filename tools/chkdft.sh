#!/bin/sh

# init
if [ $# -ne 1 -a $# -ne 2 ]; then
	echo "Usage: $0 C|FORTRAN [SINGLE|EXTENDED]"
	exit 1
fi
if [ "$2" = SINGLE ]; then
	FFTW_SFX=f
elif [ "$2" = EXTENDED ]; then
	FFTW_SFX=l
fi

# build minfft
echo Building minfft...
cc $OFLAGS -std=c99 -pedantic -Wall -Wextra -DMINFFT_$2 -c ../minfft.c
echo done.

# tests to run
dims="D1 D2 D3"
xforms="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"

# C
if [ $1 = C ]; then
	for d in $dims; do
		for x in $xforms; do
			cc chkdft.c \
			-D$d -D$x \
			-std=c99 -pedantic -Wall -Wextra \
			-DMINFFT_$2 -I .. \
			./minfft.o \
			-DFFTW_SFX=$FFTW_SFX \
			-I ~/build/fftw-mi/include \
			-L ~/build/fftw-mi/lib \
			-lfftw3$FFTW_SFX \
			-lm
			if [ $? -ne 0 ]; then
				exit
			fi
			echo "# $d $x"
			./a.out
		done
	done
fi

# Fortran
if [ $1 = FORTRAN ]; then
	gfortran -fsyntax-only -DMINFFT_$2 ../minfft.F03
	gfortran -fsyntax-only -I ~/build/fftw-mi/include fftw.f03
	for d in $dims; do
		for x in $xforms; do
			cpp chkdft.F95 -o chkdft.f95 \
			-D$d -D$x \
			-DFFTW_SFX=$FFTW_SFX \
			-P \
			&& \
			gfortran chkdft.f95 \
			-std=f2003 -pedantic -Wall -Wextra \
			./minfft.o \
			-L ~/build/fftw-mi/lib \
			-lfftw3$FFTW_SFX
			if [ $? -ne 0 ]; then
				exit
			fi
			echo "# $d $x"
			./a.out
		done
	done
fi

# cleanup
rm -f ./a.out minfft.o chkdft.f95 minfft.mod fftw.mod
