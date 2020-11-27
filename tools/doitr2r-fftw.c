#include "bench-user.h"
#include <math.h>
#include <fftw3.h>
#include <string.h>

#define CONCAT(prefix, name) prefix ## name
#if defined(BENCHFFT_SINGLE)
#define FFTW(x) CONCAT(fftwf_, x)
#elif defined(BENCHFFT_LDOUBLE)
#define FFTW(x) CONCAT(fftwl_, x)
#else
#define FFTW(x) CONCAT(fftw_, x)
#endif

BEGIN_BENCH_DOC
BENCH_DOC("name", "fftw3-r2r")
END_BENCH_DOC 

FFTW(plan) the_plan=0;
unsigned the_flags=FFTW_ESTIMATE;
fftw_r2r_kind the_kind=-1;

void useropt(const char *arg)
{
	if (!strcmp(arg,"dct2"))
		the_kind=FFTW_REDFT10;
	else if (!strcmp(arg,"dst2"))
		the_kind=FFTW_RODFT10;
	else if (!strcmp(arg,"dct3"))
		the_kind=FFTW_REDFT01;
	else if (!strcmp(arg,"dst3"))
		the_kind=FFTW_RODFT01;
	else if (!strcmp(arg,"dct4"))
		the_kind=FFTW_REDFT11;
	else if (!strcmp(arg,"dst4"))
		the_kind=FFTW_RODFT11;
}

int can_do(struct problem *p)
{
	return (p->kind==PROBLEM_REAL);
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
	int i;
	fftw_r2r_kind kind[MAX_RANK];
	for (i=0; i<p->rank; ++i)
		kind[i]=the_kind;
	the_plan=FFTW(plan_r2r)(p->rank, p->n, p->in, p->out, kind, the_flags);
	BENCH_ASSERT(the_plan);
}

void doit(int iter, struct problem *p)
{
	UNUSED(p);
	int i;
	for (i=0; i<iter; ++i) 
		FFTW(execute)(the_plan);
}

void done(struct problem *p)
{
	UNUSED(p);
	FFTW(destroy_plan)(the_plan);
}
