CC = gcc
CFLAGS = -std=c99 -pedantic -Wall -Ofast

libminfft.a: minfft.o
	ar rcs libminfft.a minfft.o

minfft.o: minfft.c minfft.h
	$(CC) $(CFLAGS) -c minfft.c -o minfft.o

clean:
	rm libminfft.a minfft.o
