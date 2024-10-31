module TestTridiag
    ! Fortran standard library.
    use stdlib_math, only: is_close, all_close
    use stdlib_linalg_constants, only: dp, ilp
    use stdlib_linalg, only: solve
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
                    new_unittest("Diagonal contructors", test_diagonal_constructors),  &
                    new_unittest("Diagonal linear solver", test_diagonal_solve)  &
                    ]
        return
    end subroutine collect_diagonal_testsuite

    subroutine test_diagonal_constructors(error)
        type(error_type), allocatable, intent(out) :: error
        type(Diagonal) :: A
        real(dp) :: dv(n), d
        integer(ilp) :: i

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

    subroutine test_diagonal_solve(error)
        type(error_type), allocatable, intent(out) :: error
        type(Diagonal) :: A
        real(dp) :: dv(n)
        
        ! Initialize matrix.
        call random_number(dv); A = Diagonal(dv)

        ! Solve with a single right-hand side vector.
        block
        real(dp) :: x(n), x_stdlib(n), b(n)
        ! Random rhs.
        call random_number(b)
        ! Solve with SpecialMatrices.
        x = solve(A, b)
        ! Solve with stdlib (dense).
        x_stdlib = solve(dense(A), b)
        ! Check error.
        call check(error, all_close(x, x_stdlib), &
        "Diagonal solve with a single rhs failed.")
        if (allocated(error)) return
        end block

        ! Solve with multiple right-hand side vectors.
        block
        real(dp) :: x(n, n), x_stdlib(n, n), b(n, n)
        ! Random rhs.
        call random_number(b)
        ! Solve with SpecialMatrices.
        x = solve(A, b)
        ! Solve with stdlib (dense).
        x_stdlib = solve(dense(A), b)
        ! Check error.
        call check(error, all_close(x, x_stdlib), &
        "Diagonal solve with multiple rhs failed.")
        end block
    end subroutine test_diagonal_solve

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
