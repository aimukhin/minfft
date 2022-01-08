set title "Speed of one-dimensional double precision in-place forward complex DFT of size N"
set xlabel "log_2(N)"
set ylabel "Mflops"
set xtics 1
set grid
log2(x)=log(x)/log(2)
set out "speed.svg"
set term svg
plot [] [0:] \
	"minfft.speed" using (log2($3)):4 with lines title "minfft" lw 2, \
	"kissfft.speed" using (log2($3)):4 with lines title "KissFFT", \
	"pocketfft.speed" using (log2($3)):4 with lines title "PocketFFT", \
	"fftw.speed" using (log2($3)):4 with lines title "FFTW"
unset out
unset term
