set title "Accuracy of one-dimensional single precision out-of-place forward complex DFT of size N"
set xlabel "log_2(N)"
set ylabel "Relative L_2 error"
set xtics 1
set format y "%.1e"
set grid
set key bottom right
log2(x)=log(x)/log(2)
set out "accuracy.svg"
set term svg
plot [] [0:] \
	"minfft-neon.accuracy" using (log2($3)):5 \
		with lines title "minfft neon" lc "blue", \
	"fftw3-neon.accuracy" using (log2($3)):5 \
		with lines title "FFTW neon" lc "green", \
	"ne10.accuracy" using (log2($3)):5 \
		with lines title "NE10" lc "red"
unset out
unset term
