.POSIX:

CC = gcc
CFLAGS = -std=c99 -pedantic -Wall -Wextra -Ofast
AS = as
ASFLAGS =

double:
	sed "s/REAL/double/" minfft.h.template > minfft.h
	sed "s/REAL/C_DOUBLE/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -c minfft.c -o minfft.o
	ar rcs libminfft.a minfft.o

single:
	sed "s/REAL/float/" minfft.h.template > minfft.h
	sed "s/REAL/C_FLOAT/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -c minfft.c -o minfft.o
	ar rcs libminfft.a minfft.o

sse3-double:
	sed "s/REAL/double/" minfft.h.template > minfft.h
	sed "s/REAL/C_DOUBLE/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -DEXT_DFT -c minfft.c -o minfft.o
	$(AS) $(ASFLAGS) machdep/sse3-double/rs_dft_1d.s -o rs_dft_1d.o
	$(AS) $(ASFLAGS) machdep/sse3-double/rs_invdft_1d.s -o rs_invdft_1d.o
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o

sse3-single:
	sed "s/REAL/float/" minfft.h.template > minfft.h
	sed "s/REAL/C_FLOAT/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -DEXT_DFT -DEXT_MKAUX -c minfft.c -o minfft.o
	$(AS) $(ASFLAGS) machdep/sse3-single/rs_dft_1d.s -o rs_dft_1d.o
	$(AS) $(ASFLAGS) machdep/sse3-single/rs_invdft_1d.s -o rs_invdft_1d.o
	$(CC) $(CFLAGS) -I. -c machdep/sse3-single/mkaux_dft_1d.c -o mkaux_dft_1d.o
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o mkaux_dft_1d.o

neon-single:
	sed "s/REAL/float/" minfft.h.template > minfft.h
	sed "s/REAL/C_FLOAT/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -DEXT_DFT -DEXT_MKAUX -c minfft.c -o minfft.o
	$(AS) $(ASFLAGS) machdep/neon-single/rs_dft_1d.s -o rs_dft_1d.o
	$(AS) $(ASFLAGS) machdep/neon-single/rs_invdft_1d.s -o rs_invdft_1d.o
	$(CC) $(CFLAGS) -I. -c machdep/neon-single/mkaux_dft_1d.c -o mkaux_dft_1d.o
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o mkaux_dft_1d.o

clean:
	rm -f *.o *.a
