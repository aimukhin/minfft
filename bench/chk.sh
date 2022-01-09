#!/bin/sh

# init
if [ -z $CC ]; then CC=cc; fi
CFLAGS="$CFLAGS -std=c99 -pedantic -Wall -Wextra -I.. -I$FFTW/include"
LDFLAGS=-L$FFTW/lib
LIBS=-lm
if [ "$1" = PERF ]; then
	DEFS="$DEFS -DPERF"
	shift
elif [ "$1" = ACC ]; then
	shift
else
	echo "Usage: $0 ACC|PERF [SINGLE|EXTENDED]"
	exit 1
fi
if [ "$1" = SINGLE ]; then
	DEFS="$DEFS -DMINFFT_SINGLE -DFFTW_SFX=f -DTHR=1.0E-06"
	LIBS="$LIBS -lfftw3f"
elif [ "$1" = EXTENDED ]; then
	DEFS="$DEFS -DMINFFT_EXTENDED -DFFTW_SFX=l -DTHR=1.0E-17"
	LIBS="$LIBS -lfftw3l"
else
	DEFS="$DEFS -DTHR=1.0E-14"
	LIBS="$LIBS -lfftw3"
fi

# build minfft
$CC $DEFS $CFLAGS -c ../minfft.c || exit 1

# run tests
XFORMS="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
SIZES="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2x1 4x2 8x4 16x8 32x16 64x32 128x64 256x128 512x256 1024x512 4x2x1 8x4x2 16x8x4 32x16x8 64x32x16 128x64x32 8x4x2x1 16x8x4x2 32x16x8x4 64x32x16x8"
for x in $XFORMS; do
	$CC -D$x $DEFS $CFLAGS chk.c minfft.o $LDFLAGS $LIBS || break
	for s in $SIZES; do
		out=$(./a.out $s) || { echo $x $s $out FAILED; break 2; }
		echo $x $s $out
	done
done

# cleanup
rm -f a.out minfft.o
