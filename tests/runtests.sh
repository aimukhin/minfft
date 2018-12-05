#!/bin/sh

# Compiler flags
cflags="../minfft.o -std=c99 -pedantic -Wall -lm"
cflags_fftw=$(echo -L ~/build/fftw/lib/ -lfftw3)
cflags_kiss=$(echo -I ~/build/kiss_fft/ ~/build/kiss_fft/kiss_fft*.o)

# Compile tests, run them and save results
do_tests() {
	for t in $1; do
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
	done
}

# Tests of the library itself
if false; then
tests=(INV)
dims=(D1 D2 D3)
xforms=(DFT REALDFT DCT2 DCT3 DCT4 DST4)
do_tests "${tests[*]}" "${dims[*]}" "${xforms[*]}"
tests=(PERF)
dims=(D1 D2 D3)
xforms=(DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4)
do_tests "${tests[*]}" "${dims[*]}" "${xforms[*]}"
fi

# Comparisons with FFTW
if false; then
tests=(CMP_FFTW PERF_FFTW)
dims=(D1 D2 D3)
xforms=(DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4)
do_tests "${tests[*]}" "${dims[*]}" "${xforms[*]}" "$cflags_fftw"
fi

# Comparisons with Kiss FFT
if false; then
tests=(CMP_KISS PERF_KISS)
dims=(D1 D2 D3)
xforms=(DFT INVDFT REALDFT INVREALDFT)
do_tests "${tests[*]}" "${dims[*]}" "${xforms[*]}" "$cflags_kiss"
fi
