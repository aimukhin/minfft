! A minimalistic FFT library
! Copyright (c) 2016-2022 Alexander Mukhin
! SPDX-License-Identifier: MIT

#if MINFFT_SINGLE
#define MINFFT_R C_FLOAT
#define MINFFT_C C_FLOAT_COMPLEX
#elif MINFFT_EXTENDED
#define MINFFT_R C_LONG_DOUBLE
#define MINFFT_C C_LONG_DOUBLE_COMPLEX
#else
#define MINFFT_R C_DOUBLE
#define MINFFT_C C_DOUBLE_COMPLEX
#endif

module minfft
    use,intrinsic :: iso_c_binding, only: minfft_real=>MINFFT_R, &
                                          minfft_cmpl=>MINFFT_C, &
                                          minfft_aux=>C_PTR, &
                                          C_INT,C_NULL_PTR
    implicit none
    interface
        subroutine minfft_dft(x,y,a) bind(C)
            import
            complex(minfft_cmpl),dimension(*),intent(in) :: x
            complex(minfft_cmpl),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_invdft(x,y,a) bind(C)
            import
            complex(minfft_cmpl),dimension(*),intent(in) :: x
            complex(minfft_cmpl),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_realdft(x,z,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            complex(minfft_cmpl),dimension(*),intent(out) :: z
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_invrealdft(z,y,a) bind(C)
            import
            complex(minfft_cmpl),dimension(*),intent(in) :: z
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_dct2(x,y,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_dst2(x,y,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_dct3(x,y,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_dst3(x,y,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_dct4(x,y,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        subroutine minfft_dst4(x,y,a) bind(C)
            import
            real(minfft_real),dimension(*),intent(in) :: x
            real(minfft_real),dimension(*),intent(out) :: y
            type(minfft_aux),value :: a
        end subroutine
        function minfft_mkaux_dft_1d(n) bind(C)
            import
            integer(C_INT),value :: n
            type(minfft_aux) :: minfft_mkaux_dft_1d
        end function
        function minfft_mkaux_dft_2d(n1,n2) bind(C)
            import
            integer(C_INT),value :: n1,n2
            type(minfft_aux) :: minfft_mkaux_dft_2d
        end function
        function minfft_mkaux_dft_3d(n1,n2,n3) bind(C)
            import
            integer(C_INT),value :: n1,n2,n3
            type(minfft_aux) :: minfft_mkaux_dft_3d
        end function
        function minfft_mkaux_dft(d,ns) bind(C)
            import
            integer(C_INT),value :: d
            integer(C_INT),dimension(*),intent(in) :: ns
            type(minfft_aux) :: minfft_mkaux_dft
        end function
        function minfft_mkaux_realdft_1d(n) bind(C)
            import
            integer(C_INT),value :: n
            type(minfft_aux) :: minfft_mkaux_realdft_1d
        end function
        function minfft_mkaux_realdft_2d(n1,n2) bind(C)
            import
            integer(C_INT),value :: n1,n2
            type(minfft_aux) :: minfft_mkaux_realdft_2d
        end function
        function minfft_mkaux_realdft_3d(n1,n2,n3) bind(C)
            import
            integer(C_INT),value :: n1,n2,n3
            type(minfft_aux) :: minfft_mkaux_realdft_3d
        end function
        function minfft_mkaux_realdft(d,ns) bind(C)
            import
            integer(C_INT),value :: d
            integer(C_INT),dimension(*),intent(in) :: ns
            type(minfft_aux) :: minfft_mkaux_realdft
        end function
        function minfft_mkaux_t2t3_1d(n) bind(C)
            import
            integer(C_INT),value :: n
            type(minfft_aux) :: minfft_mkaux_t2t3_1d
        end function
        function minfft_mkaux_t2t3_2d(n1,n2) bind(C)
            import
            integer(C_INT),value :: n1,n2
            type(minfft_aux) :: minfft_mkaux_t2t3_2d
        end function
        function minfft_mkaux_t2t3_3d(n1,n2,n3) bind(C)
            import
            integer(C_INT),value :: n1,n2,n3
            type(minfft_aux) :: minfft_mkaux_t2t3_3d
        end function
        function minfft_mkaux_t2t3(d,ns) bind(C)
            import
            integer(C_INT),value :: d
            integer(C_INT),dimension(*),intent(in) :: ns
            type(minfft_aux) :: minfft_mkaux_t2t3
        end function
        function minfft_mkaux_t4_1d(n) bind(C)
            import
            integer(C_INT),value :: n
            type(minfft_aux) :: minfft_mkaux_t4_1d
        end function
        function minfft_mkaux_t4_2d(n1,n2) bind(C)
            import
            integer(C_INT),value :: n1,n2
            type(minfft_aux) :: minfft_mkaux_t4_2d
        end function
        function minfft_mkaux_t4_3d(n1,n2,n3) bind(C)
            import
            integer(C_INT),value :: n1,n2,n3
            type(minfft_aux) :: minfft_mkaux_t4_3d
        end function
        function minfft_mkaux_t4(d,ns) bind(C)
            import
            integer(C_INT),value :: d
            integer(C_INT),dimension(*),intent(in) :: ns
            type(minfft_aux) :: minfft_mkaux_t4
        end function
        subroutine minfft_free_aux(a) bind(C)
            import
            type(minfft_aux),value :: a
        end subroutine
    end interface
end module
