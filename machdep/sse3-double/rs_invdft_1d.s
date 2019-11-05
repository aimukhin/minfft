# Double precision split-radix (2/4) DIF inverse FFT
# Arch: x86-64 SSE3
# ABI: SysV AMD64
#
# Part of the minfft library
# Copyright (c) 2018 Alexander Mukhin
# MIT License

# void
# rs_invdft_1d (
#		int N,			= %rdi
#		double complex *x,	= %rsi
#		double complex *t,	= %rdx
#		double complex *y,	= %rcx
#		int sy,			= %r8
#		const double complex *e	= %r9
# )

	.text
	.globl		rs_invdft_1d

rs_invdft_1d:
	# branches for the terminal cases
	cmp		$1,%rdi
	je		T1
	cmp		$2,%rdi
	je		T2
	cmp		$4,%rdi
	je		T4
	cmp		$8,%rdi
	je		T8
	# prepare constants
	movapd		sr(%rip),%xmm15		# sign flip mask for Re
	movapd		sri(%rip),%xmm14	# sign flip mask for Re and Im
	mov		%rdi,%rax
	shr		$2,%rax			# loop counter rax=N/4
	mov		%rdi,%r10
	shr		$1,%r10			# r10=N/2
	mov		%rdi,%r11
	add		%r10,%r11		# r11=3N/2
	# save pointers
	push		%rdx			# t
	push		%r9			# e
L:	# loop N/4 times (N>=16) with 2 iterations unrolled
	# next iteration uses registers +8 and marked with -/-
	#
	# compute temporary variables
	# t0 and t2
	movapd		(%rsi),%xmm0		# xmm0=x[n]
	addpd		(%rsi,%rdi,8),%xmm0	# xmm0=t0
	movapd		16(%rsi),%xmm8		# -/-
	addpd		16(%rsi,%rdi,8),%xmm8	# -/-
	movapd		(%rsi),%xmm2		# xmm2=x[n]
	subpd		(%rsi,%rdi,8),%xmm2	# xmm2=t2
	movapd		16(%rsi),%xmm10		# -/-
	subpd		16(%rsi,%rdi,8),%xmm10	# -/-
	# t1 and t3
	movapd		(%rsi,%r10,8),%xmm1	# xmm1=x[N/4+n]
	addpd		(%rsi,%r11,8),%xmm1	# xmm1=t1
	movapd		16(%rsi,%r10,8),%xmm9	# -/-
	addpd		16(%rsi,%r11,8),%xmm9	# -/-
	movapd		(%rsi,%r10,8),%xmm3	# xmm3=x[N/4+n]
	subpd		(%rsi,%r11,8),%xmm3	# xmm3=x[N/4+n]-x[3N/4+n]
	shufpd		$1,%xmm3,%xmm3		# swap Re and Im
	xorpd		%xmm15,%xmm3		# xmm3=t3
	movapd		16(%rsi,%r10,8),%xmm11	# -/-
	subpd		16(%rsi,%r11,8),%xmm11	# -/-
	shufpd		$1,%xmm11,%xmm11	# -/-
	xorpd		%xmm15,%xmm11		# -/-
	# store t[n]
	movapd		%xmm0,(%rdx)		# t[n]=t0
	movapd		%xmm8,16(%rdx)		# -/-
	# prepare t2+t3 and t2-t3
	movapd		%xmm2,%xmm4		# xmm4=t2
	addpd		%xmm3,%xmm2		# xmm2=t2+t3
	subpd		%xmm3,%xmm4		# xmm4=t2-t3
	movapd		%xmm10,%xmm12		# -/-
	addpd		%xmm11,%xmm10		# -/-
	subpd		%xmm11,%xmm12		# -/-
	# store t[N/4+n]
	movapd		%xmm1,(%rdx,%r10,8)	# t[N/4+n]=t1
	movapd		%xmm9,16(%rdx,%r10,8)	# -/-
	# compute t[N/2+n] and t[3N/4+n]
	#
	# 2 complex multiplications:
	# e[2n]=a,b in xmm0,1 by xmm2=c,d; res=xmm0
	# e[2n+1]=p,q in xmm3,5 by xmm4=r,s; res=xmm3
	# and 2 more multiplications for the unrolled iteration
	movddup		(%r9),%xmm0		# xmm0=a,a
	movddup		8(%r9),%xmm1		# xmm1=b,b
	mulpd		%xmm2,%xmm0		# xmm0=ac,ad
	movddup		32(%r9),%xmm8		# -/-
	movddup		40(%r9),%xmm9		# -/-
	mulpd		%xmm10,%xmm8		# -/-
	movddup		16(%r9),%xmm3		# xmm3=p,p
	movddup		24(%r9),%xmm5		# xmm5=q,q
	mulpd		%xmm4,%xmm3		# xmm3=pr,ps
	movddup		48(%r9),%xmm11		# -/-
	movddup		56(%r9),%xmm13		# -/-
	mulpd		%xmm12,%xmm11		# -/-
	shufpd		$1,%xmm2,%xmm2		# xmm2=d,c
	mulpd		%xmm2,%xmm1		# xmm1=bd,bc
	xorpd		%xmm14,%xmm1		# change sign of xmm1
	addsubpd	%xmm1,%xmm0		# xmm0=ac-bd,ad+bc
	shufpd		$1,%xmm10,%xmm10	# -/-
	mulpd		%xmm10,%xmm9		# -/-
	xorpd		%xmm14,%xmm9		# -/-
	addsubpd	%xmm9,%xmm8		# -/-
	# store t[N/2+n]
	movapd		%xmm0,(%rdx,%rdi,8)	# t[N/2+n]=xmm0
	movapd		%xmm8,16(%rdx,%rdi,8)	# -/-
	shufpd		$1,%xmm4,%xmm4		# xmm4=s,r
	mulpd		%xmm4,%xmm5		# xmm5=qs,qr
	xorpd		%xmm14,%xmm5		# change sign of xmm5
	addsubpd	%xmm5,%xmm3		# xmm3=pr-qs,ps+qr
	shufpd		$1,%xmm12,%xmm12	# -/-
	mulpd		%xmm12,%xmm13		# -/-
	xorpd		%xmm14,%xmm13		# -/-
	addsubpd	%xmm13,%xmm11		# -/-
	# store t[3N/4+n]
	movapd		%xmm3,(%rdx,%r11,8)	# t[3N/4+n]=xmm3
	movapd		%xmm11,16(%rdx,%r11,8)	# -/-
	# advance pointers
	add		$32,%rsi		# x
	add		$32,%rdx		# t
	add		$64,%r9			# e
	# decrement and check loop counter
	sub		$2,%rax
	jnz		L
	# restore pointers
	pop		%r9			# e
	pop		%rdx			# t
	# prepare arguments for the first recursive call
	lea		(%r9,%rdi,8),%r9	# r9=e+N/2
	shr		$1,%rdi			# rdi=N/2
	shl		$1,%r8			# r8=2sy
	mov		%rdx,%rsi		# x=t
	push		%rdi
	push		%rdx
	push		%rcx
	push		%r8
	push		%r9
	# call rs_invdft_1d(N/2,t,t,y,2sy,e+N/2)
	call		rs_invdft_1d
	pop		%r9			# r9=e+N/2
	pop		%r8			# r8=2sy
	pop		%rcx			# rcx=y
	pop		%rdx			# rdx=t
	pop		%rdi			# rdi=N/2
	# prepare arguments for the second recursive call
	lea		(%r9,%rdi,8),%r9	# r9=e+N/2+N/4=e+3N/4
	lea		(%rcx,%r8,8),%rcx	# rcx=y+sy
	shl		$1,%r8			# r8=4sy
	lea		(%rdx,%rdi,8),%rdx	# rdx=t+N/4
	lea		(%rdx,%rdi,8),%rdx	# rdx=t+N/4+N/4=t+N/2
	mov		%rdx,%rsi		# x=t
	shr		$1,%rdi			# rdi=N/4
	push		%rdi
	push		%rdx
	push		%rcx
	push		%r8
	push		%r9
	# call rs_invdft_1d(N/4,t+N/2,t+N/2,y+sy,4sy,e+3N/4)
	call		rs_invdft_1d
	pop		%r9			# r9=e+3N/4
	pop		%r8			# r8=4sy
	pop		%rcx			# rcx=y+sy
	pop		%rdx			# rdx=t+N/2
	pop		%rdi			# rdi=N/4
	# prepare arguments for the third recursive call
	lea		(%rcx,%r8,8),%rcx	# rcx=y+sy+2sy=y+3sy
	lea		(%rdx,%rdi,8),%rdx	# rdx=t+N/2+N/8
	lea		(%rdx,%rdi,8),%rdx	# rdx=t+N/2+N/8+N/8=t+3N/4
	mov		%rdx,%rsi		# x=t
	# no need to save the registers, since this is the last call
	# call rs_invdft_1d(N/4,t+3N/4,t+3N/4,y+3sy,4sy,e+3N/4)
	jmp		rs_invdft_1d
	# we're done

	# terminal cases
T1:	# N==1
	movapd		(%rsi),%xmm0		# xmm0=x[0]
	movapd		%xmm0,(%rcx)		# store y[0]
	ret

T2:	# N==2
	shl		$4,%r8			# r8=16sy
	movapd		(%rsi),%xmm0		# xmm0=x[0]
	movapd		16(%rsi),%xmm1		# xmm1=x[1]
	movapd		%xmm0,%xmm2
	addpd		%xmm1,%xmm0		# xmm0=x[0]+x[1]
	subpd		%xmm1,%xmm2		# xmm2=x[0]-x[1]
	movapd		%xmm0,(%rcx)		# store y[0]
	movapd		%xmm2,(%rcx,%r8)	# store y[sy]
	ret

T4:	# N==4
	# compute temporary values t0,t1,t2,t3
	movapd		(%rsi),%xmm0		# xmm0=x[0]
	movapd		32(%rsi),%xmm1		# xmm1=x[2]
	movapd		%xmm0,%xmm2		# xmm2=x[0]
	addpd		%xmm1,%xmm0		# xmm0=t0
	subpd		%xmm1,%xmm2		# xmm2=t2
	movapd		16(%rsi),%xmm1		# xmm1=x[1]
	movapd		48(%rsi),%xmm5		# xmm5=x[3]
	movapd		%xmm1,%xmm3		# xmm3=x[1]
	addpd		%xmm5,%xmm1		# xmm1=t1
	subpd		%xmm5,%xmm3		# xmm3=x[1]-x[3]
	shufpd		$1,%xmm3,%xmm3		# swap Re and Im
	xorpd		sr(%rip),%xmm3		# xmm3=t3
	# compute and store outputs
	shl		$4,%r8			# r8=16sy
	movapd		%xmm0,%xmm4		# xmm4=t0
	movapd		%xmm2,%xmm5		# xmm5=t2
	addpd		%xmm1,%xmm0		# xmm0=y[0]
	addpd		%xmm3,%xmm2		# xmm2=y[sy]
	subpd		%xmm1,%xmm4		# xmm4=y[2*sy]
	subpd		%xmm3,%xmm5		# xmm5=y[3*sy]
	movapd		%xmm0,(%rcx)		# store y[0]
	movapd		%xmm2,(%rcx,%r8)	# store y[sy]
	lea		(%rcx,%r8,2),%rcx	# advance pointer
	movapd		%xmm4,(%rcx)		# store y[2*sy]
	movapd		%xmm5,(%rcx,%r8)	# store y[3*sy]
	ret

T8:	# N==8
	# compute temporary values t00,t01,t02,t03
	movapd		(%rsi),%xmm0		# xmm0=x[0]
	movapd		64(%rsi),%xmm1		# xmm1=x[4]
	movapd		%xmm0,%xmm2		# xmm2=x[0]
	addpd		%xmm1,%xmm0		# xmm0=t0
	subpd		%xmm1,%xmm2		# xmm2=t2
	movapd		32(%rsi),%xmm1		# xmm1=x[2]
	movapd		96(%rsi),%xmm5		# xmm5=x[6]
	movapd		%xmm1,%xmm3		# xmm3=x[2]
	addpd		%xmm5,%xmm1		# xmm1=t1
	subpd		%xmm5,%xmm3		# xmm3=x[2]-x[6]
	shufpd		$1,%xmm3,%xmm3		# swap Re and Im
	xorpd		sr(%rip),%xmm3		# xmm3=t3
	movapd		%xmm0,%xmm4		# xmm4=t0
	movapd		%xmm2,%xmm5		# xmm5=t2
	addpd		%xmm1,%xmm0		# xmm0=t00
	addpd		%xmm3,%xmm2		# xmm2=t01
	subpd		%xmm1,%xmm4		# xmm4=t02
	subpd		%xmm3,%xmm5		# xmm5=t03
	# compute temporary values t10,t11,t12,t13
	movapd		16(%rsi),%xmm8		# xmm8=x[1]
	movapd		80(%rsi),%xmm9		# xmm9=x[5]
	movapd		%xmm8,%xmm10		# xmm10=x[1]
	addpd		%xmm9,%xmm8		# xmm8=t0
	subpd		%xmm9,%xmm10		# xmm10=t2
	movapd		48(%rsi),%xmm9		# xmm9=x[3]
	movapd		112(%rsi),%xmm13	# xmm13=x[7]
	movapd		%xmm9,%xmm11		# xmm11=x[3]
	addpd		%xmm13,%xmm9		# xmm9=t1
	subpd		%xmm13,%xmm11		# xmm11=x[2]-x[6]
	shufpd		$1,%xmm11,%xmm11	# swap Re and Im
	xorpd		sr(%rip),%xmm11		# xmm11=t3
	movapd		%xmm8,%xmm12		# xmm12=t0
	movapd		%xmm10,%xmm13		# xmm13=t2
	addpd		%xmm9,%xmm8		# xmm8=t10
	addpd		%xmm11,%xmm10		# xmm10=t11
	subpd		%xmm9,%xmm12		# xmm12=t12
	subpd		%xmm11,%xmm13		# xmm13=t13
	# multiply %xmm10 (t11) by e1, result in %xmm9, temp %xmm1
	movddup		e1(%rip),%xmm9		# xmm9=a,a
	movddup		e1+8(%rip),%xmm1	# xmm1=b,b
	mulpd		%xmm10,%xmm9		# xmm9=ac,ad
	shufpd		$1,%xmm10,%xmm10	# xmm10=d,c
	mulpd		%xmm10,%xmm1		# xmm1=bd,bc
	addsubpd	%xmm1,%xmm9		# xmm9=ac-bd,ad+bc
	# multiply %xmm12 (t12) by I
	shufpd		$1,%xmm12,%xmm12	# swap Re and Im
	xorpd		sr(%rip),%xmm12		# change sign of Re
	# multiply %xmm13 (t13) by e3, result in %xmm14, temp %xmm1
	movddup		e3(%rip),%xmm14		# xmm14=a,a
	movddup		e3+8(%rip),%xmm1	# xmm1=b,b
	mulpd		%xmm13,%xmm14		# xmm14=ac,ad
	shufpd		$1,%xmm13,%xmm13	# xmm13=d,c
	mulpd		%xmm13,%xmm1		# xmm1=bd,bc
	addsubpd	%xmm1,%xmm14		# xmm14=ac-bd,ad+bc
	# copy
	movapd		%xmm0,%xmm1
	movapd		%xmm2,%xmm3
	movapd		%xmm4,%xmm6
	movapd		%xmm5,%xmm7
	# compute sums
	addpd		%xmm8,%xmm0		# xmm0=t00+t10
	addpd		%xmm9,%xmm2		# xmm2=t01+t11
	addpd		%xmm12,%xmm4		# xmm4=t02+t12
	addpd		%xmm14,%xmm5		# xmm5=t03+t13
	# store outputs
	shl		$4,%r8			# r8=16*sy
	movapd		%xmm0,(%rcx)		# store y[0]
	movapd		%xmm2,(%rcx,%r8)	# store y[sy]
	lea		(%rcx,%r8,2),%rcx	# advance pointer
	movapd		%xmm4,(%rcx)		# store y[2*sy]
	movapd		%xmm5,(%rcx,%r8)	# store y[3*sy]
	lea		(%rcx,%r8,2),%rcx	# advance pointer
	# compute differences
	subpd		%xmm8,%xmm1		# xmm1=t00-t10
	subpd		%xmm9,%xmm3		# xmm3=t01-t11
	subpd		%xmm12,%xmm6		# xmm6=t02-t12
	subpd		%xmm14,%xmm7		# xmm7=t03-t13
	# store
	movapd		%xmm1,(%rcx)		# store y[4*sy]
	movapd		%xmm3,(%rcx,%r8)	# store y[5*sy]
	lea		(%rcx,%r8,2),%rcx	# advance pointer
	movapd		%xmm6,(%rcx)		# store y[6*sy]
	movapd		%xmm7,(%rcx,%r8)	# store y[7*sy]
	ret

# constants
	.align		16
# sign-flipping mask for real part
sr:	.quad		0x8000000000000000
	.quad		0x0000000000000000
# sign-flipping mask for real and imaginary parts
sri:	.quad		0x8000000000000000
	.quad		0x8000000000000000
# E_8^1
e1:	.double		0.707106781186547524400844362104849039
	.double		0.707106781186547524400844362104849039
# E_8^3
e3:	.double		-0.707106781186547524400844362104849039
	.double		0.707106781186547524400844362104849039
