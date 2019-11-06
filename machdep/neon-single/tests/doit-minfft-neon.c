#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "minfft-neon-single")
BENCH_DOC("author", "Alexander Mukhin")
BENCH_DOC("year", "2018")
BENCH_DOC("language", "C")
BENCH_DOC("email", "alexander.i.mukhin@gmail.com")
BENCH_DOC("url", "https://github.com/aimukhin/minfft/")
END_BENCH_DOC

#include "minfft.h"

minfft_aux *a; // aux data

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

int can_do(struct problem *p)
{
	int i;
	for (i=0; i<p->rank; ++i)
		if (!power_of_two(p->n[i]))
			return 0;
	return 1;
}

void setup(struct problem *p)
{
	BENCH_ASSERT(can_do(p));
	if (p->kind==PROBLEM_COMPLEX)
		a = minfft_mkaux_dft(p->rank,p->n);
	else
		a = minfft_mkaux_realdft(p->rank,p->n);
}

void doit(int iter, struct problem *p)
{
	int i;
	if (p->kind==PROBLEM_COMPLEX)
		// complex DFT
		if (p->sign<0)
			// forward transform
			for (i=0; i<iter; ++i)
				minfft_dft(p->in,p->out,a);
		else
			// inverse transform
			for (i=0; i<iter; ++i)
				minfft_invdft(p->in,p->out,a);
	else
		// real DFT
		if (p->sign<0)
			// forward transform
			for (i=0; i<iter; ++i)
				minfft_realdft(p->in,p->out,a);
		else
			// inverse transform
			for (i=0; i<iter; ++i)
				minfft_invrealdft(p->in,p->out,a);
}

void done(struct problem *p)
{
	UNUSED(p);
	minfft_free_aux(a);
}
