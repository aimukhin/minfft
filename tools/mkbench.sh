#!/bin/sh

if [ $# -ne 3 ]; then
	echo "Usage: $0 <program> <transform> <test>"
	exit 1
fi

program=$1
xform=$2
SIZES="2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 524288 1048576"

name=$($program --info name)
precision=$($program --print-precision | cut -c 1)
for size in $SIZES
do
	problem=${xform}${size}
	case $3 in
		a) out=$($program --accuracy $problem);;
		s) out=$($program --report-benchmark --speed $problem);;
		v) out=$($program --verify $problem);;
		*) exit 2;;
	esac
	echo $name ${precision}${xform} $size $out
done
