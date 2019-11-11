set title "Accuracy of one-dimensional double precision out-of-place forward complex DFT of size N"
set xlabel "log_2(N)"
set ylabel "Relative L_2 error"
set xtics 1
set grid
set key bottom right
log2(x)=log(x)/log(2)
set out "accuracy.svg"
set term svg
plot [] [0:] \
	"minfft.accuracy" using (log2($3)):5 \
		with lines title "minfft mi" lc "blue", \
	"fftw.accuracy" using (log2($3)):5 \
		with lines title "FFTW mi" lc "green", \
	"kissfft.accuracy" using (log2($3)):5 \
		with lines title "Kiss FFT" lc "red"
unset out
unset term
