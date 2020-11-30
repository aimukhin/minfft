#!/bin/sh

# init
if [ "$1" = SINGLE ]; then
	DEFS="-DMINFFT_SINGLE -DFFTW_SFX=f -DTHR=1.0E-06"
	LIBS="-lfftw3f"
	shift
elif [ "$1" = EXTENDED ]; then
	DEFS="-DMINFFT_EXTENDED -DFFTW_SFX=l -DTHR=1.0E-17"
	LIBS="-lfftw3l"
	shift
else
	DEFS="-DTHR=1.0E-14"
	LIBS="-lfftw3"
fi
if [ "$1" = PERF ]; then
	DEFS+=" -DPERF"
fi
CFLAGS+=" -std=c99 -pedantic -Wall -Wextra"
CFLAGS+=" $(echo -I.. -I ~/build/fftw/include/ -L ~/build/fftw/lib/)"
LIBS+=" -lm"

# build minfft
cc $DEFS $CFLAGS -c ../minfft.c

# run tests
XFORMS="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
SIZES="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576 2x1 4x2 8x4 16x8 32x16 64x32 128x64 256x128 512x256 1024x512 2048x1024 4x2x1 8x4x2 16x8x4 32x16x8 64x32x16 128x64x32 256x128x64"
for x in $XFORMS; do
	cc -D$x $DEFS $CFLAGS chkperf.c ./minfft.o $LIBS
	if [ $? -ne 0 ]; then
		exit 1
	fi
	for s in $SIZES; do
		echo -n "$x $s "
		./a.out $s
		if [ $? -ne 0 ]; then
			echo Failed
			exit 1
		fi
	done
done

# cleanup
rm ./a.out minfft.o
