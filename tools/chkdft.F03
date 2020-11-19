! Compare minfft with FFTW
! Dimensionality: D1 D2 D3
! Transforms: DFT INVDFT REALDFT INVREALDFT DCT2 DST2 DCT3 DST3 DCT4 DST4

! macros for calling variants of FFTW
#ifndef FFTW_SFX
#define FFTW_SFX
#endif
#define FFTW(F) FFTW2(F,FFTW_SFX)
#define FFTW2(F,S) FFTW22(F,S)
#define FFTW22(F,S) fftw##S##_##F

program main
    use minfft
    use fftw
    implicit none
#if D1
    ! one-dimensional transforms
    integer,parameter :: maxn=65536*16 ! max transform size
    integer :: n ! current transform size
    real(minfft_real) :: r ! relative error
    type(minfft_aux) :: a ! minfft aux
    type(C_PTR) :: p ! FFTW plan
#if DFT || INVDFT
    complex(minfft_cmpl),dimension(:),allocatable :: mi,mo,fi,fo
    real(minfft_real) :: re,im
    integer :: i
#endif
#if REALDFT
    real(minfft_real),dimension(:),allocatable :: mi,fi
    complex(minfft_cmpl),dimension(:),allocatable :: mo,fo
#endif
#if INVREALDFT
    complex(minfft_cmpl),dimension(:),allocatable :: mi,fi
    real(minfft_real),dimension(:),allocatable :: mo,fo
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
    real(minfft_real),dimension(:),allocatable :: mi,mo,fi,fo
#endif
    n=1
    do while (n<=maxn)
        ! allocate and initialize arrays
#if DFT || INVDFT
        allocate(mi(n),mo(n),fi(n),fo(n))
        do i=1,n
            call random_number(re)
            call random_number(im)
            mi(i)=CMPLX(re-0.5,im-0.5,kind(mi))
        end do
        fi=mi
#endif
#if REALDFT
        allocate(mi(n),mo(n/2+1),fi(n),fo(n/2+1))
        call random_number(mi)
        mi=mi-0.5
        fi=mi
#endif
#if INVREALDFT
        allocate(mi(n/2+1),mo(n),fi(n/2+1),fo(n))
        call random_number(mo)
        mo=mo-0.5
        a=minfft_mkaux_realdft_1d(n)
        call minfft_realdft(mo,mi,a)
        fi=mi
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
        allocate(mi(n),mo(n),fi(n),fo(n))
        call random_number(mi)
        mi=mi-0.5
        fi=mi
#endif
        ! do transforms
#if DFT
        a=minfft_mkaux_dft_1d(n)
        call minfft_dft(mi,mo,a)
        p=FFTW(plan_dft_1d)(n,fi,fo,FFTW_FORWARD,FFTW_ESTIMATE)
        call FFTW(execute_dft)(p,fi,fo)
#endif
#if INVDFT
        a=minfft_mkaux_dft_1d(n)
        call minfft_invdft(mi,mo,a)
        p=FFTW(plan_dft_1d)(n,fi,fo,FFTW_BACKWARD,FFTW_ESTIMATE)
        call FFTW(execute_dft)(p,fi,fo)
#endif
#if REALDFT
        a=minfft_mkaux_realdft_1d(n)
        call minfft_realdft(mi,mo,a)
        p=FFTW(plan_dft_r2c_1d)(n,fi,fo,FFTW_ESTIMATE)
        call FFTW(execute_dft_r2c)(p,fi,fo)
#endif
#if INVREALDFT
        call minfft_invrealdft(mi,mo,a)
        p=FFTW(plan_dft_c2r_1d)(n,fi,fo,FFTW_ESTIMATE)
        call FFTW(execute_dft_c2r)(p,fi,fo)
#endif
#if DCT2
        a=minfft_mkaux_t2t3_1d(n)
        call minfft_dct2(mi,mo,a)
        p=FFTW(plan_r2r_1d)(n,fi,fo,FFTW_REDFT10,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST2
        a=minfft_mkaux_t2t3_1d(n)
        call minfft_dst2(mi,mo,a)
        p=FFTW(plan_r2r_1d)(n,fi,fo,FFTW_RODFT10,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DCT3
        a=minfft_mkaux_t2t3_1d(n)
        call minfft_dct3(mi,mo,a)
        p=FFTW(plan_r2r_1d)(n,fi,fo,FFTW_REDFT01,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST3
        a=minfft_mkaux_t2t3_1d(n)
        call minfft_dst3(mi,mo,a)
        p=FFTW(plan_r2r_1d)(n,fi,fo,FFTW_RODFT01,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DCT4
        a=minfft_mkaux_t4_1d(n)
        call minfft_dct4(mi,mo,a)
        p=FFTW(plan_r2r_1d)(n,fi,fo,FFTW_REDFT11,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST4
        a=minfft_mkaux_t4_1d(n)
        call minfft_dst4(mi,mo,a)
        p=FFTW(plan_r2r_1d)(n,fi,fo,FFTW_RODFT11,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
        ! compare results
        r=maxval(abs(mo-fo))/maxval(abs(fo))
        print *,n,r
        if (r>=T) stop 1 ! error threshold exceeded
        ! free plans
        call minfft_free_aux(a)
        call FFTW(destroy_plan)(p)
        ! free arrays
        deallocate(mi,mo,fi,fo)
        ! next transform size
        n=n*2
    end do
#endif
#if D2
    ! two-dimensional transforms
    integer,parameter :: maxn=1024 ! max transform size
    integer :: n ! current transform size
    real(minfft_real) :: r ! relative error
    type(minfft_aux) :: a ! minfft aux
    type(C_PTR) :: p ! FFTW plan
#if DFT || INVDFT
    complex(minfft_cmpl),dimension(:,:),allocatable :: mi,mo,fi,fo
    real(minfft_real) :: re,im
    integer :: i,j
#endif
#if REALDFT
    real(minfft_real),dimension(:,:),allocatable :: mi,fi
    complex(minfft_cmpl),dimension(:,:),allocatable :: mo,fo
#endif
#if INVREALDFT
    complex(minfft_cmpl),dimension(:,:),allocatable :: mi,fi
    real(minfft_real),dimension(:,:),allocatable :: mo,fo
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
    real(minfft_real),dimension(:,:),allocatable :: mi,mo,fi,fo
#endif
    n=1
    do while (n<=maxn)
        ! allocate and initialize arrays
#if DFT || INVDFT
        allocate(mi(n,2*n),mo(n,2*n),fi(n,2*n),fo(n,2*n))
        do j=1,2*n
            do i=1,n
                call random_number(re)
                call random_number(im)
                mi(i,j)=CMPLX(re-0.5,im-0.5,kind(mi))
            end do
        end do
        fi=mi
#endif
#if REALDFT
        allocate(mi(n,2*n),mo(n/2+1,2*n),fi(n,2*n),fo(n/2+1,2*n))
        call random_number(mi)
        mi=mi-0.5
        fi=mi
#endif
#if INVREALDFT
        allocate(mi(n/2+1,2*n),mo(n,2*n),fi(n/2+1,2*n),fo(n,2*n))
        call random_number(mo)
        mo=mo-0.5
        a=minfft_mkaux_realdft_2d(2*n,n)
        call minfft_realdft(mo,mi,a)
        fi=mi
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
        allocate(mi(n,2*n),mo(n,2*n),fi(n,2*n),fo(n,2*n))
        call random_number(mi)
        mi=mi-0.5
        fi=mi
#endif
        ! do transforms
#if DFT
        a=minfft_mkaux_dft_2d(2*n,n)
        call minfft_dft(mi,mo,a)
        p=FFTW(plan_dft_2d)(2*n,n,fi,fo,FFTW_FORWARD,FFTW_ESTIMATE)
        call FFTW(execute_dft)(p,fi,fo)
#endif
#if INVDFT
        a=minfft_mkaux_dft_2d(2*n,n)
        call minfft_invdft(mi,mo,a)
        p=FFTW(plan_dft_2d)(2*n,n,fi,fo,FFTW_BACKWARD,FFTW_ESTIMATE)
        call FFTW(execute_dft)(p,fi,fo)
#endif
#if REALDFT
        a=minfft_mkaux_realdft_2d(2*n,n)
        call minfft_realdft(mi,mo,a)
        p=FFTW(plan_dft_r2c_2d)(2*n,n,fi,fo,FFTW_ESTIMATE)
        call FFTW(execute_dft_r2c)(p,fi,fo)
#endif
#if INVREALDFT
        call minfft_invrealdft(mi,mo,a)
        p=FFTW(plan_dft_c2r_2d)(2*n,n,fi,fo,FFTW_ESTIMATE)
        call FFTW(execute_dft_c2r)(p,fi,fo)
#endif
#if DCT2
        a=minfft_mkaux_t2t3_2d(2*n,n)
        call minfft_dct2(mi,mo,a)
        p=FFTW(plan_r2r_2d)(2*n,n,fi,fo,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST2
        a=minfft_mkaux_t2t3_2d(2*n,n)
        call minfft_dst2(mi,mo,a)
        p=FFTW(plan_r2r_2d)(2*n,n,fi,fo,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DCT3
        a=minfft_mkaux_t2t3_2d(2*n,n)
        call minfft_dct3(mi,mo,a)
        p=FFTW(plan_r2r_2d)(2*n,n,fi,fo,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST3
        a=minfft_mkaux_t2t3_2d(2*n,n)
        call minfft_dst3(mi,mo,a)
        p=FFTW(plan_r2r_2d)(2*n,n,fi,fo,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DCT4
        a=minfft_mkaux_t4_2d(2*n,n)
        call minfft_dct4(mi,mo,a)
        p=FFTW(plan_r2r_2d)(2*n,n,fi,fo,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST4
        a=minfft_mkaux_t4_2d(2*n,n)
        call minfft_dst4(mi,mo,a)
        p=FFTW(plan_r2r_2d)(2*n,n,fi,fo,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
        ! compare results
        r=maxval(abs(mo-fo))/maxval(abs(fo))
        print *,2*n,n,r
        if (r>=T) stop 1 ! error threshold exceeded
        ! free plans
        call minfft_free_aux(a)
        call FFTW(destroy_plan)(p)
        ! free arrays
        deallocate(mi,mo,fi,fo)
        ! next transform size
        n=n*2
    end do
#endif
#if D3
    ! three-dimensional transforms
    integer,parameter :: maxn=64 ! max transform size
    integer :: n ! current transform size
    real(minfft_real) :: r ! relative error
    type(minfft_aux) :: a ! minfft aux
    type(C_PTR) :: p ! FFTW plan
#if DFT || INVDFT
    complex(minfft_cmpl),dimension(:,:,:),allocatable :: mi,mo,fi,fo
    real(minfft_real) :: re,im
    integer :: i,j,k
#endif
#if REALDFT
    real(minfft_real),dimension(:,:,:),allocatable :: mi,fi
    complex(minfft_cmpl),dimension(:,:,:),allocatable :: mo,fo
#endif
#if INVREALDFT
    complex(minfft_cmpl),dimension(:,:,:),allocatable :: mi,fi
    real(minfft_real),dimension(:,:,:),allocatable :: mo,fo
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
    real(minfft_real),dimension(:,:,:),allocatable :: mi,mo,fi,fo
#endif
    n=1
    do while (n<=maxn)
        ! allocate and initialize arrays
#if DFT || INVDFT
        allocate(mi(n,2*n,4*n),mo(n,2*n,4*n),fi(n,2*n,4*n),fo(n,2*n,4*n))
        do k=1,4*n
            do j=1,2*n
                do i=1,n
                    call random_number(re)
                    call random_number(im)
                    mi(i,j,k)=CMPLX(re-0.5,im-0.5,kind(mi))
                end do
            end do
        end do
        fi=mi
#endif
#if REALDFT
        allocate(mi(n,2*n,4*n),mo(n/2+1,2*n,4*n),fi(n,2*n,4*n),fo(n/2+1,2*n,4*n))
        call random_number(mi)
        mi=mi-0.5
        fi=mi
#endif
#if INVREALDFT
        allocate(mi(n/2+1,2*n,4*n),mo(n,2*n,4*n),fi(n/2+1,2*n,4*n),fo(n,2*n,4*n))
        call random_number(mo)
        mo=mo-0.5
        a=minfft_mkaux_realdft_3d(4*n,2*n,n)
        call minfft_realdft(mo,mi,a)
        fi=mi
#endif
#if DCT2 || DST2 || DCT3 || DST3 || DCT4 || DST4
        allocate(mi(n,2*n,4*n),mo(n,2*n,4*n),fi(n,2*n,4*n),fo(n,2*n,4*n))
        call random_number(mi)
        mi=mi-0.5
        fi=mi
#endif
        ! do transforms
#if DFT
        a=minfft_mkaux_dft_3d(4*n,2*n,n)
        call minfft_dft(mi,mo,a)
        p=FFTW(plan_dft_3d)(4*n,2*n,n,fi,fo,FFTW_FORWARD,FFTW_ESTIMATE)
        call FFTW(execute_dft)(p,fi,fo)
#endif
#if INVDFT
        a=minfft_mkaux_dft_3d(4*n,2*n,n)
        call minfft_invdft(mi,mo,a)
        p=FFTW(plan_dft_3d)(4*n,2*n,n,fi,fo,FFTW_BACKWARD,FFTW_ESTIMATE)
        call FFTW(execute_dft)(p,fi,fo)
#endif
#if REALDFT
        a=minfft_mkaux_realdft_3d(4*n,2*n,n)
        call minfft_realdft(mi,mo,a)
        p=FFTW(plan_dft_r2c_3d)(4*n,2*n,n,fi,fo,FFTW_ESTIMATE)
        call FFTW(execute_dft_r2c)(p,fi,fo)
#endif
#if INVREALDFT
        call minfft_invrealdft(mi,mo,a)
        p=FFTW(plan_dft_c2r_3d)(4*n,2*n,n,fi,fo,FFTW_ESTIMATE)
        call FFTW(execute_dft_c2r)(p,fi,fo)
#endif
#if DCT2
        a=minfft_mkaux_t2t3_3d(4*n,2*n,n)
        call minfft_dct2(mi,mo,a)
        p=FFTW(plan_r2r_3d)(4*n,2*n,n,fi,fo,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST2
        a=minfft_mkaux_t2t3_3d(4*n,2*n,n)
        call minfft_dst2(mi,mo,a)
        p=FFTW(plan_r2r_3d)(4*n,2*n,n,fi,fo,FFTW_RODFT10,FFTW_RODFT10,FFTW_RODFT10,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DCT3
        a=minfft_mkaux_t2t3_3d(4*n,2*n,n)
        call minfft_dct3(mi,mo,a)
        p=FFTW(plan_r2r_3d)(4*n,2*n,n,fi,fo,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST3
        a=minfft_mkaux_t2t3_3d(4*n,2*n,n)
        call minfft_dst3(mi,mo,a)
        p=FFTW(plan_r2r_3d)(4*n,2*n,n,fi,fo,FFTW_RODFT01,FFTW_RODFT01,FFTW_RODFT01,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DCT4
        a=minfft_mkaux_t4_3d(4*n,2*n,n)
        call minfft_dct4(mi,mo,a)
        p=FFTW(plan_r2r_3d)(4*n,2*n,n,fi,fo,FFTW_REDFT11,FFTW_REDFT11,FFTW_REDFT11,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
#if DST4
        a=minfft_mkaux_t4_3d(4*n,2*n,n)
        call minfft_dst4(mi,mo,a)
        p=FFTW(plan_r2r_3d)(4*n,2*n,n,fi,fo,FFTW_RODFT11,FFTW_RODFT11,FFTW_RODFT11,FFTW_ESTIMATE)
        call FFTW(execute_r2r)(p,fi,fo)
#endif
        ! compare results
        r=maxval(abs(mo-fo))/maxval(abs(fo))
        print *,4*n,2*n,n,r
        if (r>=T) stop 1 ! error threshold exceeded
        ! free plans
        call minfft_free_aux(a)
        call FFTW(destroy_plan)(p)
        ! free arrays
        deallocate(mi,mo,fi,fo)
        ! next transform size
        n=n*2
    end do
#endif
end program
