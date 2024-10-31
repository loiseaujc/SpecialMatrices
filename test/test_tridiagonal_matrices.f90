module TestTridiag
    ! Fortran standard library.
    use stdlib_kinds, only: dp
    use stdlib_math, only: is_close, all_close
    ! Testdrive.
    use testdrive, only: new_unittest, unittest_type, error_type, check
    ! SpecialMatrices
    use SpecialMatrices
    implicit none
    private

    integer, parameter :: n = 512

    public :: collect_diagonal_testsuite

contains

    !-------------------------------------
    !-----     DIAGONAL MATRICES     -----
    !-------------------------------------

    subroutine collect_diagonal_testsuite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
                    new_unittest("Diagonal contructors", test_diagonal_constructors) &
                    ]
        return
    end subroutine collect_diagonal_testsuite

    subroutine test_diagonal_constructors(error)
        type(error_type), allocatable, intent(out) :: error
        type(Diagonal), allocatable :: A
        real(dp) :: dv(n), d
        integer :: i

        ! Zero-initialization of Diagonal.
        A = Diagonal(n)
        call check(error, all(A%dv == [(0.0_dp, i=1, n)]), &
        "Initializing the zero Diagonal matrix failed.")
        if (allocated(error)) return

        ! Array-based initialization of Diagonal.
        dv = [(i, i=1, n)] ; A = Diagonal(dv)
        call check(error, all(A%dv == dv), &
        "Array-based initialization of Diagonal matrix failed.")
        if (allocated(error)) return

        ! Constant-Diagonal matrix.
        d = 1.0_dp ; A = Diagonal(d, n)
        call check(error, all(A%dv == [(d, i=1, n)]), &
        "Constant-Diagonal initialization failed.")

        return
    end subroutine test_diagonal_constructors

    !---------------------------------------
    !-----     BIDIAGONAL MATRICES     -----
    !---------------------------------------

    !----------------------------------------
    !-----     TRIDIAGONAL MATRICES     -----
    !----------------------------------------

    !--------------------------------------------------
    !-----     SYMMETRIC TRIDIAGONAL MATRICES     -----
    !--------------------------------------------------

end module
