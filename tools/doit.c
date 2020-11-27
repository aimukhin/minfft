#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "minfft")
BENCH_DOC("version", "1.2.0")
BENCH_DOC("year", "2020")
BENCH_DOC("author", "Alexander Mukhin")
BENCH_DOC("language", "C")
BENCH_DOC("email", "alexander.i.mukhin@gmail.com")
BENCH_DOC("url", "https://github.com/aimukhin/minfft/")
BENCH_DOC("copyright",
"MIT License\n"
"\n"
"Copyright (c) 2020 Alexander Mukhin\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining a copy\n"
"of this software and associated documentation files (the \"Software\"), to deal\n"
"in the Software without restriction, including without limitation the rights\n"
"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
"copies of the Software, and to permit persons to whom the Software is\n"
"furnished to do so, subject to the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be included in all\n"
"copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n"
"SOFTWARE.\n"
)
END_BENCH_DOC

#include "minfft.h"

minfft_aux *a; // aux data

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p,out,-1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p,in,-1.0);
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
		a=minfft_mkaux_dft(p->rank,p->n);
	else
		a=minfft_mkaux_realdft(p->rank,p->n);
	BENCH_ASSERT(a!=NULL);
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
