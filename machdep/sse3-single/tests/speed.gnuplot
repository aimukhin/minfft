set title "Speed of one-dimensional single precision out-of-place forward complex DFT of size N"
set xlabel "log_2(N)"
set ylabel "Mflops"
set xtics 1
set ytics 500
set grid
log2(x)=log(x)/log(2)
set out "speed.svg"
set term svg
plot [] [0:] \
	"minfft-sse3.speed" using (log2($3)):4 \
		with lines title "minfft sse3-single" lc "blue", \
	"fftw3-sse3.speed" using (log2($3)):4 \
		with lines title "FFTW sse2/sse3" lc "green"
unset out
unset term
