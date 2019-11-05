@ Single precision split-radix (2/4) DIF inverse FFT
@ Arch: ARM NEON
@ ABI: ARM Architecture Procedure Call Standard (AAPCS)
@
@ Part of the minfft library
@ Copyright (c) 2018 Alexander Mukhin
@ MIT License

@ void
@ rs_invdft_1d (
@		int N,			= r0
@		float complex *x,	= r1
@		float complex *t,	= r2
@		float complex *y,	= r3
@		int sy,			= [sp]
@		const float complex *e	= [sp,#4]
@ )

	.text
	.fpu		neon
	.globl		rs_invdft_1d

rs_invdft_1d:
	@ select branches
	cmp		r0,#1
	beq		T1
	cmp		r0,#2
	beq		T2
	cmp		r0,#4
	beq		T4
	cmp		r0,#8
	beq		T8
	@ prepare for loop
	push		{r4-r11,lr}		@ save registers
	lsr		r4,r0,#2		@ r4 counts from N/4 to 1
	ldr		r5,[sp,#40]		@ r5=e
	add		r6,r1,r0,lsl#1		@ r6=x+N/4
	add		r7,r1,r0,lsl#2		@ r7=x+N/2
	add		r8,r6,r0,lsl#2		@ r8=x+3N/4
	add		r9,r2,r0,lsl#1		@ r9=t+N/4
	add		r10,r2,r0,lsl#2		@ r10=t+N/2
	add		r11,r9,r0,lsl#2		@ r11=t+3N/4
L:	@ main loop
	@ two iterations unrolled, ' means next iteration value
	@ load x
	vld2.32		{q0},[r1]!		@ q0=x[n]_r,x[n']_r,x[n]_i,x[n']_i
	vld2.32		{q1},[r6]!		@ q1=..[n+n/4]..
	vld2.32		{q2},[r7]!		@ q2=..[n+N/2]..
	vld2.32		{q3},[r8]!		@ q3=..[n+3N/4]..
	@ compute temporary variables t0,t1,t2,t3
	vadd.f32	q8,q0,q2		@ q8=t0_r,t0'_r,t0_i,t0'_i
	vadd.f32	q9,q1,q3		@ q9=..t1..
	vsub.f32	q10,q0,q2		@ q10=..t2..
	vsub.f32	d22,d7,d3
	vsub.f32	d23,d2,d6		@ q11=..t3..
	@ store t[n] and t[n+N/4]
	vst2.32		{q8},[r2]!		@ store t[n] and t[n']
	vst2.32		{q9},[r9]!		@ ..[n+N/4]..
	@ prepare t2+t3 and t2-t3
	vadd.f32	d0,d20,d22		@ d0=(t2+t3)_r,(t2+t3)'_r
	vsub.f32	d1,d20,d22		@ d1=(t2-t3)_r,(t2-t3)'_r
	vadd.f32	d2,d21,d23		@ d2=(t2+t3)_i,(t2+t3)'_i
	vsub.f32	d3,d21,d23		@ d3=(t2-t3)_i,(t2-t3)'_i
	@ load and conjugate exponents
	vld1.32		{q2-q3},[r5]!
	vneg.f32	q3,q3
	@ compute products: (t2+t3)*e1, (t2+t3)'*e1', (t2-t3)*e3, (t2-t3)'*e3'
	vmul.f32	q8,q0,q2
	vmul.f32	q9,q0,q3
	vmls.f32	q8,q1,q3
	vmla.f32	q9,q1,q2
	@ store t[n+N/2] and t[n+3N/4]
	vst2.32		{d16,d18},[r10]!	@ store t[n+N/2] and '
	vst2.32		{d17,d19},[r11]!	@ ..[n+3N/4]..
	@ decrement loop counter
	subs		r4,r4,#2
	@ branch if not zero
	bne		L
	@ restore pointers
	sub		r2,r2,r0,lsl#1		@ t
	sub		r5,r5,r0,lsl#2		@ e
	@ fetch parameter
	ldr		r4,[sp,#36]		@ r4=sy
	@ recursive calls
	lsr		r0,r0,#1		@ N/2
	mov		r1,r2			@ x=t
	lsl		r4,r4,#1		@ 2sy
	add		r5,r5,r0,lsl#3		@ e+N/2
	push		{r0-r3}			@ save registers
	push		{r4-r5}			@ pass sy and e
	@ call rs_invdft_1d(N/2,t,t,y,2*sy,e+N/2)
	bl		rs_invdft_1d
	pop		{r4-r5}			@ restore sy and e
	pop		{r0-r3}			@ restore other parameters
	lsr		r0,r0,#1		@ N/4
	lsl		r4,r4,#1		@ 4sy
	add		r2,r2,r0,lsl#4		@ t+N/2
	mov		r1,r2			@ x=t
	add		r3,r3,r4,lsl#1		@ y+sy
	add		r5,r5,r0,lsl#3		@ e+N/2+N/4
	push		{r0-r3}			@ save registers
	push		{r4-r5}			@ pass sy and e
	@ call rs_invdft_1d(N/4,t+N/2,t+N/2,y+sy,4*sy,e+3*N/4)
	bl		rs_invdft_1d
	pop		{r4-r5}			@ restore sy and e
	pop		{r0-r3}			@ restore other parameters
	add		r2,r2,r0,lsl#3		@ t+N/2+N/4
	mov		r1,r2			@ x=t
	add		r3,r3,r4,lsl#2		@ y+sy+2sy
	push		{r4-r5}			@ pass sy and e
	@ call rs_invdft_1d(N/4,t+3*N/4,t+3*N/4,y+3*sy,4*sy,e+3*N/4)
	bl		rs_invdft_1d
	add		sp,sp,#8		@ drop parameters
	pop		{r4-r11,lr}		@ restore ARM registers
	@ done
	bx		lr

	@ terminal branches
T1:	vld1.32		d0,[r1]
	vst1.32		d0,[r3]
	bx		lr

T2:	@ load data
	vld1.32		{d0,d1},[r1]
	@ prepare offset
	ldr		r2,[sp]			@ r2=sy
	lsl		r2,r2,#3
	@ add, subtract, and store
	vadd.f32	d2,d0,d1
	vst1.32		d2,[r3],r2		@ store y[0]
	vsub.f32	d3,d0,d1
	vst1.32		d3,[r3]			@ store y[sy]
	@ done
	bx		lr

T4:	@ load data
	vld1.32		{d0,d1,d2,d3},[r1]
	@ compute temporary variables
	vadd.f32	q2,q0,q1		@ d4=t0, d5=t1
	vsub.f32	d6,d0,d2		@ d6=t2
	vsub.f32	s14,s7,s3
	vsub.f32	s15,s2,s6		@ d7=t3
	@ prepare offset
	ldr		r2,[sp]			@ r2=sy
	lsl		r2,r2,#3
	@ add, subtract, and store
	vadd.f32	d0,d4,d5
	vst1.32		d0,[r3],r2		@ y[0]=t0+t1
	vadd.f32	d1,d6,d7
	vst1.32		d1,[r3],r2		@ y[sy]=t2+t3
	vsub.f32	d2,d4,d5
	vst1.32		d2,[r3],r2		@ y[2sy]=t0-t1
	vsub.f32	d3,d6,d7
	vst1.32		d3,[r3]			@ y[3sy]=t2-t3
	@ done
	bx		lr

T8:	@ save NEON registers
	vpush		{q4-q7}
	@ load data
	vld1.32		{d0,d1,d2,d3},[r1]!
	vld1.32		{d4,d5,d6,d7},[r1]
	@ compute temporary variables
	vadd.f32	q4,q0,q2		@ q4=t0,t0'
	vadd.f32	q5,q1,q3		@ q5=t1,t1'
	vsub.f32	q6,q0,q2		@ q6=t2,t2'
	vsub.f32	s29,s4,s12
	vsub.f32	s28,s13,s5
	vsub.f32	s31,s6,s14
	vsub.f32	s30,s15,s7		@ q7=t3,t3'
	vadd.f32	q0,q4,q5		@ q0=t00,t10
	vadd.f32	q1,q6,q7		@ q1=t01,t11
	vsub.f32	d4,d8,d10
	vsub.f32	s10,s23,s19
	vsub.f32	s11,s18,s22		@ q2=t02,I*t12
	vsub.f32	q3,q6,q7		@ q3=t03,t13
	@ multiply t11 by e1
	ldr		r2,=e
	vld1.32		{d8},[r2]!		@ d8=e1
	vmul.f32	s18,s6,s16
	vmul.f32	s19,s6,s17
	vmls.f32	s18,s7,s17
	vmla.f32	s19,s7,s16
	vmov.f32	d3,d9
	@ multiply t13 by e3
	vld1.32		{d8},[r2]		@ d8=e3
	vmul.f32	s18,s14,s16
	vmul.f32	s19,s14,s17
	vmls.f32	s18,s15,s17
	vmla.f32	s19,s15,s16
	vmov.f32	d7,d9
	@ prepare offset
	ldr		r2,[sp,#64]		@ r2=sy
	lsl		r2,r2,#3
	@ add and store
	vadd.f32	d8,d0,d1
	vst1.32		d8,[r3],r2		@ store y[0]
	vadd.f32	d9,d2,d3
	vst1.32		d9,[r3],r2		@ store y[sy]
	vadd.f32	d10,d4,d5
	vst1.32		d10,[r3],r2		@ store y[2*sy]
	vadd.f32	d11,d6,d7
	vst1.32		d11,[r3],r2		@ store y[3*sy]
	@ subtract and store
	vsub.f32	d8,d0,d1
	vst1.32		d8,[r3],r2		@ store y[4*sy]
	vsub.f32	d9,d2,d3
	vst1.32		d9,[r3],r2		@ store y[5*sy]
	vsub.f32	d10,d4,d5
	vst1.32		d10,[r3],r2		@ store y[6*sy]
	vsub.f32	d11,d6,d7
	vst1.32		d11,[r3]		@ store y[7*sy]
	@ restore NEON registers
	vpop		{q4-q7}
	@ done
	bx		lr

@ constants
	.data
e:	.float		0.707106781186547524400844362104849039
	.float		0.707106781186547524400844362104849039
	.float		-0.707106781186547524400844362104849039
	.float		0.707106781186547524400844362104849039
