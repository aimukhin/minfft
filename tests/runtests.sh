#!/bin/sh

# Compiler flags
cflags="-std=c99 -pedantic -Wall -lm -L.. -lminfft"
cflags_fftw=$(echo -L ~/build/fftw/lib/ -lfftw3)
cflags_kiss=$(echo -I ~/build/kiss_fft/ ~/build/kiss_fft/kiss_fft*.o)

# Compile tests, run them and save results
do_tests() {
	t=$1
	for d in $2; do
		for x in $3; do
			echo -n $t $d $x
			out=results/${t}_${d}_${x}
			if [ ! -e $out ]; then
				echo
				cc -D$t -D$d -D$x chkdft.c $cflags $4
				echo "# $t $d $x" > $out
				taskset -c 1 ./a.out >> $out
			else
				echo " *** SKIP"
			fi
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
	dims=(D1 D2 D3)
	xforms=(DFT REALDFT DCT2 DCT3 DCT4 DST4)
	do_tests INV "${dims[*]}" "${xforms[*]}"
fi
## Performance
if [ $t = PERF -o $t = ALL ]; then
	dims=(D1 D2 D3)
	xforms=(DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4)
	do_tests PERF "${dims[*]}" "${xforms[*]}"
fi

# Comparisons with FFTW
## Accuracy
if [ $t = CMP_FFTW -o $t = ALL ]; then
	dims=(D1 D2 D3)
	xforms=(DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4)
	do_tests CMP_FFTW "${dims[*]}" "${xforms[*]}" "$cflags_fftw"
fi
## Performance
if [ $t = PERF_FFTW -o $t = ALL ]; then
	dims=(D1 D2 D3)
	xforms=(DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4)
	do_tests PERF_FFTW "${dims[*]}" "${xforms[*]}" "$cflags_fftw"
fi

# Comparisons with Kiss FFT
## Accuracy
if [ $t = CMP_KISS -o $t = ALL ]; then
	dims=(D1 D2 D3)
	xforms=(DFT INVDFT REALDFT INVREALDFT)
	do_tests CMP_KISS "${dims[*]}" "${xforms[*]}" "$cflags_kiss"
fi
## Performance
if [ $t = PERF_KISS -o $t = ALL ]; then
	dims=(D1 D2 D3)
	xforms=(DFT INVDFT REALDFT INVREALDFT)
	do_tests PERF_KISS "${dims[*]}" "${xforms[*]}" "$cflags_kiss"
fi
