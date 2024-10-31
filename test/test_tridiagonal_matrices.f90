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
                    new_unittest("Diagonal matmul", test_diagonal_matmul),  &
                    new_unittest("Diagonal linear solver", test_diagonal_solve),  &
                    new_unittest("Bidiagonal matmul", test_bidiagonal_matmul),  &
                    new_unittest("Bidiagonal linear solver", test_bidiagonal_solve),  &
                    new_unittest("Tridiagonal matmul", test_tridiagonal_matmul),  &
                    new_unittest("Tridiagonal linear solver", test_tridiagonal_solve),  &
                    new_unittest("SymTridiagonal matmul", test_symtridiagonal_matmul),  &
                    new_unittest("SymTridiagonal linear solver", test_symtridiagonal_solve)  &
                    ]
        return
    end subroutine collect_diagonal_testsuite

    subroutine test_diagonal_matmul(error)
        type(error_type), allocatable, intent(out) :: error
        type(Diagonal) :: A
        real(dp) :: dv(n)

        ! Initialize matrix.
        call random_number(dv) ; A = Diagonal(dv)

        ! Matrix-vector product.
        block
        real(dp) :: x(n), y(n), y_dense(n)
        call random_number(x)
        y = matmul(A, x) ; y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Diagonal matrix-vector product failed.")
        if (allocated(error)) return
        end block

        ! Matrix-matrix product.
        block
        real(dp) :: x(n, n), y(n, n), y_dense(n, n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Diagonal matrix-matrix product failed.")
        end block
        return
    end subroutine test_diagonal_matmul

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

        return
    end subroutine test_diagonal_solve

    !---------------------------------------
    !-----     BIDIAGONAL MATRICES     -----
    !---------------------------------------

    subroutine test_bidiagonal_matmul(error)
        type(error_type), allocatable, intent(out) :: error
        type(Bidiagonal) :: A
        real(dp) :: dv(n), ev(n-1)

        ! Initialize matrix.
        call random_number(dv); call random_number(ev)
        A = Bidiagonal(dv, ev)

        ! Matrix-vector product with A lower bidiagonal.
        block
        real(dp) :: x(n), y(n), y_dense(n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Lower-bidiagonal matrix-vector product failed.")
        if (allocated(error)) return

        A%which = "U"
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Upper-bidiagonal matrix-vector product failed.")
        if (allocated(error)) return
        end block

        ! Matrix-matrix product with A lower bidiagonal.
        block
        real(dp) :: x(n, n), y(n, n), y_dense(n, n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Upper-bidiagonal matrix-matrix product failed.")
        if (allocated(error)) return

        A%which = "L"
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Lower-bidiagonal matrix-matrix product failed.")
        end block
        return
    end subroutine test_Bidiagonal_matmul

    subroutine test_bidiagonal_solve(error)
        type(error_type), allocatable, intent(out) :: error
        type(Bidiagonal) :: A
        real(dp) :: ev(n-1), dv(n)

        ! Initialize matrix.
        call random_number(ev); call random_number(dv)
        A = Bidiagonal(dv, ev)

        ! Solve with a singe right-hand side vector.
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
        "Lower-bidiagonal solve with a single rhs failed.")
        if (allocated(error)) return

        A%which = "U"
        ! Solve with SpecialMatrices.
        x = solve(A, b)
        ! Solve with stdlib (dense).
        x_stdlib = solve(dense(A), b)
        ! Check error.
        call check(error, all_close(x, x_stdlib), &
        "Upper-bidiagonal solve with a single rhs failed.")
        if (allocated(error)) return
        end block

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
        "Upper-bidiagonal solve with multiple rhs failed.")
        if(allocated(error)) return

        A%which = "L"
        ! Solve with SpecialMatrices.
        x = solve(A, b)
        ! Solve with stdlib (dense).
        x_stdlib = solve(dense(A), b)
        ! Check error.
        call check(error, all_close(x, x_stdlib), &
        "Lower-bidiagonal solve with multiple rhs failed.")
        end block

        return
    end subroutine test_bidiagonal_solve

    !----------------------------------------
    !-----     TRIDIAGONAL MATRICES     -----
    !----------------------------------------

    subroutine test_tridiagonal_matmul(error)
        type(error_type), allocatable, intent(out) :: error
        type(Tridiagonal) :: A
        real(dp) :: dl(n-1), dv(n), du(n-1)

        ! Initialize matrix.
        call random_number(dl); call random_number(dv); call random_number(du)
        A = Tridiagonal(dl, dv, du)

        ! Matrix-vector product.
        block
        real(dp) :: x(n), y(n), y_dense(n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Tridiagonal matrix-vector product failed.")
        if (allocated(error)) return
        end block

        ! Matrix-matrix product.
        block
        real(dp) :: x(n, n), y(n, n), y_dense(n, n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "Tridiagonal matrix-matrix product failed.")
        end block
        return
    end subroutine test_tridiagonal_matmul

    subroutine test_tridiagonal_solve(error)
        type(error_type), allocatable, intent(out) :: error
        type(Tridiagonal) :: A
        real(dp) :: dl(n-1), dv(n), du(n-1)

        ! Initialize matrix.
        call random_number(dl); call random_number(dv); call random_number(du)
        A = Tridiagonal(dl, dv, du)

        ! Solve with a singe right-hand side vector.
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
        "Tridiagonal solve with a single rhs failed.")
        if (allocated(error)) return
        end block

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
        "Tridiagonal solve with multiple rhs failed.")
        end block

        return
    end subroutine test_tridiagonal_solve

    !--------------------------------------------------
    !-----     SYMMETRIC TRIDIAGONAL MATRICES     -----
    !--------------------------------------------------

    subroutine test_symtridiagonal_matmul(error)
        type(error_type), allocatable, intent(out) :: error
        type(SymTridiagonal) :: A
        real(dp) :: dv(n), ev(n-1)

        ! Initialize matrix.
        call random_number(dv); call random_number(ev)
        A = SymTridiagonal(dv, ev)

        ! Matrix-vector product.
        block
        real(dp) :: x(n), y(n), y_dense(n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "SymTridiagonal matrix-vector product failed.")
        if (allocated(error)) return
        end block

        ! Matrix-matrix product.
        block
        real(dp) :: x(n, n), y(n, n), y_dense(n, n)
        call random_number(x)
        y = matmul(A, x); y_dense = matmul(dense(A), x)
        call check(error, all_close(y, y_dense), &
        "SymTridiagonal matrix-matrix product failed.")
        end block
        return
    end subroutine test_symtridiagonal_matmul

    subroutine test_symtridiagonal_solve(error)
        type(error_type), allocatable, intent(out) :: error
        type(SymTridiagonal) :: A
        real(dp) :: dv(n), ev(n-1)

        ! Initialize matrix.
        call random_number(dv); call random_number(ev)
        A = SymTridiagonal(dv, ev)

        ! Solve with a singe right-hand side vector.
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
        "SymTridiagonal solve with a single rhs failed.")
        if (allocated(error)) return
        end block

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
        "SymTridiagonal solve with multiple rhs failed.")
        end block

        return
    end subroutine test_symtridiagonal_solve

end module
