#!/bin/sh

# init
if [ "$1" = SINGLE ]; then
	FFTW_SFX=f
elif [ "$1" = EXTENDED ]; then
	FFTW_SFX=l
elif [ -n "$1" ]; then
	echo "Usage: $0 [SINGLE|EXTENDED]"
	exit 1
fi

# build minfft
cc $CFLAGS -std=c99 -pedantic -Wall -Wextra -DMINFFT_$1 -c ../minfft.c

# tests to run
XFORMS="DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4"
SIZES="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576"

for x in $XFORMS
do
	# build chkperf
	cc -std=c99 -Wall -Wextra -pedantic chkperf.c -DMINFFT_$1 -DFFTW_SFX=$FFTW_SFX -I.. -D$x -lm ./minfft.o -I ~/build/fftw/include/ -L ~/build/fftw/lib/ -lfftw3$FFTW_SFX
	if [ $? -ne 0 ]; then
		exit 1
	fi
	# run tests
	for s in $SIZES
	do
		echo "$x $s $(./a.out $s)"
	done
done

# cleanup
rm ./a.out minfft.o
