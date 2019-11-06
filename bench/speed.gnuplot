set title "Speed of one-dimensional double precision out-of-place forward complex DFT of size N"
set xlabel "log_2(N)"
set ylabel "Mflops"
set xtics 1
set ytics 200
set grid
log2(x)=log(x)/log(2)
set out "speed.svg"
set term svg
plot [] [0:] \
	"minfft-mi.speed" using (log2($3)):4 \
		with lines title "minfft mi" lc "blue", \
	"fftw3-mi.speed" using (log2($3)):4 \
		with lines title "FFTW mi" lc "green", \
	"kissfft.speed" using (log2($3)):4 \
		with lines title "Kiss FFT" lc "red"
unset out
unset term
