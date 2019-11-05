set title "Speed of one-dimensional single precision out-of-place forward complex DFT of size N"
set xlabel "log_2(N)"
set ylabel "Mflops"
set xtics 1
set ytics 100
set grid
set key bottom right
log2(x)=log(x)/log(2)
set out "speed.svg"
set term svg
plot [] [0:] \
	"minfft-neon-single.speed" using (log2($3)):4 \
		with lines title "minfft (neon-single branch)" lc "blue", \
	"fftw3-neon.speed" using (log2($3)):4 \
		with lines title "FFTW (neon optimizations)" lc "green", \
	"ne10.speed" using (log2($3)):4 \
		with lines title "NE10" lc "red"
unset out
unset term
