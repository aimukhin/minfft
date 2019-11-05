#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "Ne10")
END_BENCH_DOC

#include "NE10.h"

void *cfg; // Ne10 config

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
	return (p->rank==1 && power_of_two(p->n[0]));
}

void setup(struct problem *p)
{
	BENCH_ASSERT(can_do(p));
	ne10_init();
	if (p->kind==PROBLEM_COMPLEX)
		cfg = ne10_fft_alloc_c2c_float32(p->n[0]);
	else
		cfg = ne10_fft_alloc_r2c_float32(p->n[0]);
}

void doit(int iter, struct problem *p)
{
	int i;
	if (p->kind==PROBLEM_COMPLEX)
		// complex DFT
		if (p->sign<0)
			// forward transform
			for (i=0; i<iter; ++i)
				ne10_fft_c2c_1d_float32(p->out,p->in,cfg,0);
		else
			// inverse transform
			for (i=0; i<iter; ++i)
				ne10_fft_c2c_1d_float32(p->out,p->in,cfg,1);
	else
		// real DFT
		if (p->sign<0)
			// forward transform
			for (i=0; i<iter; ++i)
				ne10_fft_r2c_1d_float32(p->out,p->in,cfg);
		else
			// inverse transform
			for (i=0; i<iter; ++i)
				ne10_fft_c2r_1d_float32(p->out,p->in,cfg);
}

void done(struct problem *p)
{
	if (p->kind==PROBLEM_COMPLEX)
		ne10_fft_destroy_c2c_float32(cfg);
	else
		ne10_fft_destroy_r2c_float32(cfg);
}
