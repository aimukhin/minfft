.POSIX:

CC = gcc
CFLAGS = -std=c99 -pedantic -Wall -Wextra -Ofast
AS = as
ASFLAGS =

double:
	sed "s/TYPE/double/" minfft.h.template > minfft.h
	sed "s/TYPE/C_DOUBLE/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -c minfft.c
	ar rcs libminfft.a minfft.o

single:
	sed "s/TYPE/float/" minfft.h.template > minfft.h
	sed "s/TYPE/C_FLOAT/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -c minfft.c
	ar rcs libminfft.a minfft.o

extended:
	sed "s/TYPE/long double/" minfft.h.template > minfft.h
	sed "s/TYPE/C_LONG_DOUBLE/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -c minfft.c
	ar rcs libminfft.a minfft.o

sse3-double:
	sed "s/TYPE/double/" minfft.h.template > minfft.h
	sed "s/TYPE/C_DOUBLE/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -DEXT_DFT -c minfft.c
	$(AS) $(ASFLAGS) machdep/sse3-double/rs_dft_1d.s -o rs_dft_1d.o
	$(AS) $(ASFLAGS) machdep/sse3-double/rs_invdft_1d.s -o rs_invdft_1d.o
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o

sse3-single:
	sed "s/TYPE/float/" minfft.h.template > minfft.h
	sed "s/TYPE/C_FLOAT/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -DEXT_DFT -DEXT_MKAUX -c minfft.c
	$(AS) $(ASFLAGS) machdep/sse3-single/rs_dft_1d.s -o rs_dft_1d.o
	$(AS) $(ASFLAGS) machdep/sse3-single/rs_invdft_1d.s -o rs_invdft_1d.o
	$(CC) $(CFLAGS) -I. -c machdep/sse3-single/mkaux_dft_1d.c
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o mkaux_dft_1d.o

neon-single:
	sed "s/TYPE/float/" minfft.h.template > minfft.h
	sed "s/TYPE/C_FLOAT/" minfft.f03.template > minfft.f03
	$(CC) $(CFLAGS) -DEXT_DFT -DEXT_MKAUX -c minfft.c
	$(AS) $(ASFLAGS) machdep/neon-single/rs_dft_1d.s -o rs_dft_1d.o
	$(AS) $(ASFLAGS) machdep/neon-single/rs_invdft_1d.s -o rs_invdft_1d.o
	$(CC) $(CFLAGS) -I. -c machdep/neon-single/mkaux_dft_1d.c
	ar rcs libminfft.a minfft.o rs_dft_1d.o rs_invdft_1d.o mkaux_dft_1d.o

clean:
	rm -f minfft.h minfft.f03 *.o *.a
