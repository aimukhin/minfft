#!/bin/sh

if [ $# -ne 2 ]; then
	echo "Usage: $0 <program> <test>"
	exit 1
fi

program=$1
reality=c
place=o
direction=f
SIZES="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576"

name=$($program --info name)
precision=$($program --print-precision | cut -c 1)
for size in $SIZES
do
	problem=${reality}${place}${direction}${size}
	case $2 in
		a) out=$($program --accuracy $problem);;
		s) out=$($program --report-benchmark --speed $problem);;
		v) out=$($program --verify $problem);;
		*) exit 2;;
	esac
	echo $name ${precision}${reality}${place}${direction} $size $out
done
