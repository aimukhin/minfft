# Single precision split-radix (2/4) DIF inverse FFT
# Arch: x86-64 SSE3
# ABI: SysV AMD64
#
# Part of the minfft library
# Copyright (c) 2019 Alexander Mukhin
# MIT License

# void
# rs_invdft_1d (
#		int N,			= %rdi
#		float complex *x,	= %rsi
# 		float complex *t,	= %rdx
#		float complex *y,	= %rcx
#		int sy,			= %r8
#		const float complex *e	= %r9
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
	movaps		lhr(%rip),%xmm15	# sign mask (lo and hi Re)
	movaps		lhri(%rip),%xmm14	# sign mask (lo and hi Re and Im)
	mov		%rdi,%rax
	shr		$2,%rax			# loop counter rax=N/4
	mov		%rdi,%r10
	shr		$1,%r10			# r10=N/2
	mov		%rdi,%r11
	add		%r10,%r11		# r11=3N/2
	# save pointers
	push		%rdx			# t
	push		%r9			# e
L:	# loop N/4 times (N>=16) with 2 iterations in parallel
	# and two iterations unrolled
	# ' means next iteration value
	# unrolled iteration uses registers +8 and marked with -/-
	#
	# compute temporary variables
	# t0,t0' and t2,t2'
	movaps		(%rsi),%xmm0		# xmm0=x[n],'
	addps		(%rsi,%rdi,4),%xmm0	# xmm0=t0,'
	movaps		(%rsi),%xmm2		# xmm2=x[n],'
	subps		(%rsi,%rdi,4),%xmm2	# xmm2=t2,'
	movaps		16(%rsi),%xmm8
	addps		16(%rsi,%rdi,4),%xmm8
	movaps		16(%rsi),%xmm10
	subps		16(%rsi,%rdi,4),%xmm10
	# t1,t1' and t3,t3'
	movaps		(%rsi,%r10,4),%xmm1	# xmm1=x[N/4+n],'
	addps		(%rsi,%r11,4),%xmm1	# xmm1=t1,'
	movaps		(%rsi,%r10,4),%xmm3	# xmm3=x[N/4+n],'
	subps		(%rsi,%r11,4),%xmm3	# xmm3=x[N/4+n]-x[3N/4+n],'
	shufps		$0xB1,%xmm3,%xmm3	# swap Re and Im
	xorps		%xmm15,%xmm3		# xmm3=t3,'
	movaps		16(%rsi,%r10,4),%xmm9
	addps		16(%rsi,%r11,4),%xmm9
	movaps		16(%rsi,%r10,4),%xmm11
	subps		16(%rsi,%r11,4),%xmm11
	shufps		$0xB1,%xmm11,%xmm11
	xorps		%xmm15,%xmm11
	# store t[n] and t[n']
	movaps		%xmm0,(%rdx)		# t[n],t[n']=xmm0
	movaps		%xmm8,16(%rdx)
	# prepare t2-t3,t2'-t3' and t2+t3,t2'+t3'
	movaps		%xmm2,%xmm4		# xmm4=t2,'
	addps		%xmm3,%xmm2		# xmm2=t2+t3,'
	subps		%xmm3,%xmm4		# xmm4=t2-t3,'
	movaps		%xmm10,%xmm12
	addps		%xmm11,%xmm10
	subps		%xmm11,%xmm12
	# store t[N/4+n] and t[N/4+n']
	movaps		%xmm1,(%rdx,%r10,4)	# t[N/4+n],t[N/4+n']=xmm1
	movaps		%xmm9,16(%rdx,%r10,4)
	# compute t[N/2+n] and t[3N/4+n]
	#
	# 2 complex multiplications (each consists of two in parallel)
	# e[2n]=er,ei,e'r,e'i in xmm0,1 by xmm2=cr,ci,c'r,c'i; res=xmm0
	# e[2n+1]=gr,gi,g'r,g'i in xmm3,5 by xmm4=dr,di,d'r,d'i; res=xmm3
	movsldup	(%r9),%xmm0		# xmm0=er,er,'
	movshdup	(%r9),%xmm1		# xmm1=ei,ei,'
	mulps		%xmm2,%xmm0		# xmm0=cr*er,ci*er,'
	movsldup	32(%r9),%xmm8
	movshdup	32(%r9),%xmm9
	mulps		%xmm10,%xmm8
	movsldup	16(%r9),%xmm3		# xmm3=gr,gr,'
	movshdup	16(%r9),%xmm5		# xmm5=gi,gi,'
	mulps		%xmm4,%xmm3		# xmm3=dr*gr,di*gr,'
	movsldup	48(%r9),%xmm11
	movshdup	48(%r9),%xmm13
	mulps		%xmm12,%xmm11
	shufps		$0xB1,%xmm2,%xmm2	# xmm2=ci,cr,'
	mulps		%xmm2,%xmm1		# xmm1=ci*ei,cr*ei,'
	xorps		%xmm14,%xmm1		# change sign
	addsubps	%xmm1,%xmm0		# xmm0=cr*er-ci*ei,ci*er+cr*ei,'
	shufps		$0xB1,%xmm10,%xmm10
	mulps		%xmm10,%xmm9
	xorps		%xmm14,%xmm9		# change sign
	addsubps	%xmm9,%xmm8
	# store t[N/2+n] and t[N/2+n']
	movaps		%xmm0,(%rdx,%rdi,4)	# t[N/2+n],t[N/2+n']=xmm0
	movaps		%xmm8,16(%rdx,%rdi,4)
	shufps		$0xB1,%xmm4,%xmm4	# xmm4=di,dr,'
	mulps		%xmm4,%xmm5		# xmm5=di*gi,dr*gi,'
	xorps		%xmm14,%xmm5		# change sign
	addsubps	%xmm5,%xmm3		# xmm3=dr*gr-di*gi,di*gr+dr*gi,'
	shufps		$0xB1,%xmm12,%xmm12
	mulps		%xmm12,%xmm13
	xorps		%xmm14,%xmm13		# change sign
	addsubps	%xmm13,%xmm11
	# store t[3N/4+n] and t[3N/4+n']
	movaps		%xmm3,(%rdx,%r11,4)	# t[3N/4+n],t[3N/4+n']=xmm3
	movaps		%xmm11,16(%rdx,%r11,4)
	# advance pointers
	add		$32,%rsi		# x
	add		$32,%rdx		# t
	add		$64,%r9			# e
	# decrement and check loop counter
	sub		$4,%rax
	jnz		L
	# restore pointers
	pop		%r9			# e
	pop		%rdx			# t
	# prepare arguments for the first recursive call
	lea		(%r9,%rdi,4),%r9	# r9=e+N/2
	shr		$1,%rdi			# rdi=N/2
	shl		$1,%r8			# r8=2sy
	mov		%rdx,%rsi		# x=t
	push		%rdi
	push		%rdx
	push		%rcx
	push		%r8
	push		%r9
	# call rs_invdft_1d(N/2,x,y,2sy,e+N/2)
	call		rs_invdft_1d
	pop		%r9			# r9=e+N/2
	pop		%r8			# r8=2sy
	pop		%rcx			# rcx=y
	pop		%rdx			# rdx=t
	pop		%rdi			# rdi=N/2
	# prepare arguments for the second recursive call
	lea		(%r9,%rdi,4),%r9	# r9=e+N/2+N/4=e+3N/4
	lea		(%rcx,%r8,4),%rcx	# rcx=y+sy
	shl		$1,%r8			# r8=4sy
	lea		(%rdx,%rdi,4),%rdx	# rdx=t+N/4
	lea		(%rdx,%rdi,4),%rdx	# rdx=t+N/4+N/4=t+N/2
	mov		%rdx,%rsi		# x=t
	shr		$1,%rdi			# rdi=N/4
	push		%rdi
	push		%rdx
	push		%rcx
	push		%r8
	push		%r9
	# call rs_invdft_1d(N/4,x+N/2,y+sy,4sy,e+3N/4)
	call		rs_invdft_1d
	pop		%r9			# r9=e+3N/4
	pop		%r8			# r8=4sy
	pop		%rcx			# rcx=y+sy
	pop		%rdx			# rdx=t+N/2
	pop		%rdi			# rdi=N/4
	# prepare arguments for the third recursive call
	lea		(%rcx,%r8,4),%rcx	# rcx=y+sy+2sy=y+3sy
	lea		(%rdx,%rdi,4),%rdx	# rdx=t+N/2+N/8
	lea		(%rdx,%rdi,4),%rdx	# rdx=t+N/2+N/8+N/8=t+3N/4
	mov		%rdx,%rsi		# x=t
	# no need to save the registers, since this is the last call
	# call rs_invdft_1d(N/4,x+3N/4,y+3sy,4sy,e+3N/4)
	jmp		rs_invdft_1d
	# we're done

	# terminal cases
T1:	# N==1
	movlps		(%rsi),%xmm0		# xmm0=x[0]
	movlps		%xmm0,(%rcx)		# store y[0]
	ret

T2:	# N==2
	movlps		(%rsi),%xmm0		# xmm0=x[0]
	movlps		8(%rsi),%xmm1		# xmm1=x[1]
	movaps		%xmm0,%xmm2
	addps		%xmm1,%xmm0		# xmm0=x[0]+x[1]
	subps		%xmm1,%xmm2		# xmm2=x[0]-x[1]
	movlps		%xmm0,(%rcx)		# store y[0]
	movlps		%xmm2,(%rcx,%r8,8)	# store y[sy]
	ret

T4:	# N==4
	# compute temporary values t0,t1,t2,t3
	movaps		(%rsi),%xmm0		# xmm0=x[0],x[1]
	addps		16(%rsi),%xmm0		# xmm0=t0,t1
	movaps		(%rsi),%xmm2		# xmm2=x[0],x[1]
	subps		16(%rsi),%xmm2		# xmm2=t2,almost-t3
	shufps		$0xB4,%xmm2,%xmm2	# swap Re and Im of high complex
	xorps		hr(%rip),%xmm2		# xmm2=t2,t3
	movaps		%xmm2,%xmm1		# xmm1=t2,t3
	movhlps		%xmm0,%xmm1		# xmm1=t1,t3
	movlhps		%xmm2,%xmm0		# xmm0=t0,t2
	movaps		%xmm0,%xmm2		# xmm2=t0,t2
	addps		%xmm1,%xmm0		# xmm0=t0+t1,t2+t3
	subps		%xmm1,%xmm2		# xmm2=t0-t1,t2-t3
	# compute and store outputs
	shl		$3,%r8			# r8=8sy
	movlps		%xmm0,(%rcx)		# y[0]=t0+t1
	movhps		%xmm0,(%rcx,%r8)	# y[sy]=t2+t3
	lea		(%rcx,%r8,2),%rcx	# advance pointer
	movlps		%xmm2,(%rcx)		# y[2*sy]=t0-t1
	movhps		%xmm2,(%rcx,%r8)	# y[3*sy]=t2-t3
	ret

T8:	# N==8
	# compute temporary variables t00,t01,t02,t03 and t10,t11,t12,t13
	movlps		(%rsi),%xmm0		# xmm0=x[0]
	movhps		16(%rsi),%xmm0		# xmm0=x[0],x[2]
	movlps		32(%rsi),%xmm1		# xmm1=x[4]
	movhps		48(%rsi),%xmm1		# xmm1=x[6]
	movaps		%xmm0,%xmm2		# xmm2=x[0],x[2]
	addps		%xmm1,%xmm0		# xmm0=t0,t1
	subps		%xmm1,%xmm2		# xmm2=t2,almost-t3
	movlps		8(%rsi),%xmm4		# xmm4=x[1]
	movhps		24(%rsi),%xmm4		# xmm4=x[1],x[3]
	movlps		40(%rsi),%xmm5		# xmm5=x[5]
	movhps		56(%rsi),%xmm5		# xmm5=x[7]
	movaps		%xmm4,%xmm6		# xmm6=x[1],x[3]
	addps		%xmm5,%xmm4		# xmm4=t0',t1'
	subps		%xmm5,%xmm6		# xmm6=t2',almost-t3'
	shufps		$0xB4,%xmm2,%xmm2	# swap Re and Im of high complex
	xorps		hr(%rip),%xmm2		# xmm2=t2,t3
	shufps		$0xB4,%xmm6,%xmm6	# swap Re and Im of high complex
	xorps		hr(%rip),%xmm6		# xmm6=t2',t3'
	movaps		%xmm2,%xmm1		# xmm1=t2,t3
	movhlps		%xmm0,%xmm1		# xmm1=t1,t3
	movlhps		%xmm2,%xmm0		# xmm0=t0,t2
	movaps		%xmm0,%xmm2		# xmm2=t0,t2
	movaps		%xmm6,%xmm5		# xmm5=t2',t3'
	movhlps		%xmm4,%xmm5		# xmm5=t1',t3'
	movlhps		%xmm6,%xmm4		# xmm4=t0',t2'
	movaps		%xmm4,%xmm6		# xmm6=t0',t2'
	addps		%xmm1,%xmm0		# xmm0=t0+t1,t2+t3=t00,t01
	subps		%xmm1,%xmm2		# xmm2=t0-t1,t2-t3=t02,t03
	addps		%xmm5,%xmm4		# xmm4=t0'+t1',t2'+t3'=t10,t11
	subps		%xmm5,%xmm6		# xmm6=t0'-t1',t2'-t3'=t12,t13
	# rearrange for complex multiplications
	movapd		%xmm6,%xmm7		# xmm7=t12,t13
	movhlps		%xmm4,%xmm7		# xmm7=t11,t13
	# multiply t11,t13 (xmm7) by e1,e3
	movsldup	e(%rip),%xmm8		# xmm8=e1r:e1r,e3r:e3r
	movshdup	e(%rip),%xmm9		# xmm9=e1i:e1i,e3i:e3i
	mulps		%xmm7,%xmm8		# xmm8=t11r*e1r:t11i*e1r,'
	shufps		$0xB1,%xmm7,%xmm7	# xmm7=t11i:t11r,t13i:t13r
	mulps		%xmm7,%xmm9		# xmm9=t11i*e1i:t11r*e1i,'
	addsubps	%xmm9,%xmm8		# xmm8=t11*e1,t13*e3
	# multiply t12 by I
	shufps		$0xE1,%xmm6,%xmm6	# swap low Re and Im
	xorps		lr(%rip),%xmm6		# change sign of low Re
	# rearrange
	movlhps		%xmm8,%xmm4		# xmm4=t10,t11*e1
	shufps		$0xE4,%xmm8,%xmm6	# xmm6=t12,t13*e3
	# add, subtract, and store
	movaps		%xmm0,%xmm1		# xmm1=t00,t03
	movaps		%xmm2,%xmm3		# xmm3=x02,t01
	addps		%xmm4,%xmm0		# xmm0=t00+t10,t01+t11
	addps		%xmm6,%xmm2		# xmm2=t02+t12,t03+t13
	subps		%xmm4,%xmm1		# xmm1=t00-t10,t01-t11
	subps		%xmm6,%xmm3		# xmm3=t02-t12,t03-t13
	movlps		%xmm0,(%rcx)		# store y[0]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movhps		%xmm0,(%rcx)		# store y[sy]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movlps		%xmm2,(%rcx)		# store y[2*sy]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movhps		%xmm2,(%rcx)		# store y[3*sy]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movlps		%xmm1,(%rcx)		# store y[4*sy]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movhps		%xmm1,(%rcx)		# store y[5*sy]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movlps		%xmm3,(%rcx)		# store y[6*sy]
	lea		(%rcx,%r8,8),%rcx	# advance pointer
	movhps		%xmm3,(%rcx)		# store y[7*sy]
	# done
	ret

# constants
	.align		16
# sign-flipping mask for low and high real parts
lhr:	.long		0x80000000
	.long		0x00000000
	.long		0x80000000
	.long		0x00000000
# sign-flipping mask for low and high real and imaginary parts
lhri:	.long		0x80000000
	.long		0x80000000
	.long		0x80000000
	.long		0x80000000
# sign-flipping mask for high real part
hr:	.long		0x00000000
	.long		0x00000000
	.long		0x80000000
	.long		0x00000000
# sign-flipping mask for low real part
lr:	.long		0x80000000
	.long		0x00000000
	.long		0x00000000
	.long		0x00000000
# E_8^1, E_8^3
e:	.float		0.707106781186547524400844362104849039
	.float		0.707106781186547524400844362104849039
	.float		-0.707106781186547524400844362104849039
	.float		0.707106781186547524400844362104849039
