#include "bench-user.h"
#include <string.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "minfft-r2r")
END_BENCH_DOC 

#include "minfft.h"

minfft_aux *a;
enum {DCT2=0,DST2,DCT3,DST3,DCT4,DST4} the_kind=-1;

void useropt(const char *arg)
{
	if (!strcmp(arg,"dct2"))
		the_kind=DCT2;
	else if (!strcmp(arg,"dst2"))
		the_kind=DST2;
	else if (!strcmp(arg,"dct3"))
		the_kind=DCT3;
	else if (!strcmp(arg,"dst3"))
		the_kind=DST3;
	else if (!strcmp(arg,"dct4"))
		the_kind=DCT4;
	else if (!strcmp(arg,"dst4"))
		the_kind=DST4;
}

int can_do(struct problem *p)
{
	int i;
	if (p->kind!=PROBLEM_REAL)
		return 0;
	for (i=0; i<p->rank; ++i)
		if (!power_of_two(p->n[i]))
			return 0;
	return 1;
}

void copy_h2c(struct problem *p, bench_complex *out)
{
	copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
	copy_c2h_1d_halfcomplex(p, in, -1.0);
}

void setup(struct problem *p)
{
	BENCH_ASSERT(can_do(p));
	BENCH_ASSERT(the_kind!=-1);
	if (the_kind==DCT2||the_kind==DST2||the_kind==DCT3||the_kind==DST3)
		a=minfft_mkaux_t2t3(p->rank,p->n);
	else
		a=minfft_mkaux_t4(p->rank,p->n);
	BENCH_ASSERT(a!=NULL);
}

void doit(int iter, struct problem *p)
{
	int i;
	if (the_kind==DCT2)
		for (i=0; i<iter; ++i)
			minfft_dct2(p->in,p->out,a);
	else if (the_kind==DST2)
		for (i=0; i<iter; ++i)
			minfft_dst2(p->in,p->out,a);
	else if (the_kind==DCT3)
		for (i=0; i<iter; ++i)
			minfft_dct3(p->in,p->out,a);
	else if (the_kind==DST3)
		for (i=0; i<iter; ++i)
			minfft_dst3(p->in,p->out,a);
	else if (the_kind==DCT4)
		for (i=0; i<iter; ++i)
			minfft_dct4(p->in,p->out,a);
	else if (the_kind==DST4)
		for (i=0; i<iter; ++i)
			minfft_dst4(p->in,p->out,a);
}

void done(struct problem *p)
{
	UNUSED(p);
	minfft_free_aux(a);
}
