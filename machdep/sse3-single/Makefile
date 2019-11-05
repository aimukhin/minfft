CC = gcc
CFLAGS = -std=c99 -pedantic -Wall -Ofast
AS = as
ASFLAGS =

libminfft.a: minfft.o rs_dft_1d.o rs_invdft_1d.o
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o

minfft.o: ../../minfft.c ../../minfft.h
	$(CC) $(CFLAGS) -DMACHDEP -I../.. -c ../../minfft.c -o minfft.o

rs_dft_1d.o: rs_dft_1d.s
	$(AS) $(ASFLAGS) rs_dft_1d.s -o rs_dft_1d.o

rs_invdft_1d.o: rs_invdft_1d.s
	$(AS) $(ASFLAGS) rs_invdft_1d.s -o rs_invdft_1d.o

clean:
	rm libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o
