module test_symtridiagonal
   ! Fortran standard library.
   use stdlib_math, only: is_close, all_close
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg, only: diag, det, trace, inv, solve, svdvals, eigvalsh
   ! Testdrive.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   ! SpecialMatrices
   use SpecialMatrices
   implicit none
   private

   integer, parameter :: n = 512

   public :: collect_symtridiagonal_testsuite

contains

   !--------------------------------------------------
   !-----     SYMMETRIC TRIDIAGONAL MATRICES     -----
   !--------------------------------------------------
   subroutine collect_symtridiagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("SymTridiagonal scalar multiplication", test_scalar_multiplication), &
                  new_unittest("SymTridiagonal trace", test_trace), &
                  new_unittest("SymTridiagonal determinant", test_det), &
                  new_unittest("SymTridiagonal matmul", test_matmul), &
                  new_unittest("SymTridiagonal linear solver", test_solve), &
                  new_unittest("Pos. def. SymTridiagonal linear solver", test_posdef_solve), &
                  new_unittest("SymTridiagonal eigenvalue decomposition", test_eigh) &
                  ]
      return
   end subroutine collect_symtridiagonal_testsuite

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1))
      call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "SymTridiagonal matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "SymTridiagonal matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1))
      call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Solve with a singe right-hand side vector.
      block
         real(dp), allocatable :: x(:), x_stdlib(:), b(:)
         allocate (b(n))
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
         real(dp), allocatable :: x(:, :), x_stdlib(:, :), b(:, :)
         allocate (b(n, n))
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
   end subroutine test_solve

   subroutine test_posdef_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A

      ! Initialize matrix.
      A = SymTridiagonal(2.0_dp, -1.0_dp, n, .true.)

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
                    "Pos. def. SymTridiagonal solve with a single rhs failed.")
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
                    "Pos. def. SymTridiagonal solve with multiple rhs failed.")
      end block

      return
   end subroutine test_posdef_solve

   subroutine test_scalar_multiplication(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A, B
      real(dp), allocatable :: dv(:), ev(:)
      real(dp) :: alpha

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1)); call random_number(dv); call random_number(ev)
      A = SymTriDiagonal(dv, ev)
      call random_number(alpha)

      ! Scalar-matrix multiplication.
      B = alpha*A
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "alpha*SymTridiagonal failed.")
      if (allocated(error)) return

      ! Matrix-scalar multipliation.
      B = A*alpha
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "SymTridiagonal*alpha failed.")
      return
   end subroutine test_scalar_multiplication

   subroutine test_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1)); call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "SymTridiagonal trace failed.")
      return
   end subroutine test_trace

   subroutine test_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1)); call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "SymTridiagonal det failed.")
      return
   end subroutine test_det

   subroutine test_eigh(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)
      real(dp), allocatable :: lambda(:), vectors(:, :), Amat(:, :)
      real(dp), allocatable :: lambda_stdlib(:)
      integer :: i, j

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1)); call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Compute eigenvalues.
      lambda = eigvalsh(A); lambda_stdlib = eigvalsh(dense(A))

      ! Check error.
      call check(error, all_close(lambda, lambda_stdlib), &
                 "SymTridiagonal eigvalsh failed.")
      if (allocated(error)) return

      ! Compute eigenvalues and eigenvectors.
      call eigh(A, lambda, vectors)

      ! Check error.
      Amat = matmul(vectors, matmul(diag(lambda), transpose(vectors)))
      call check(error, maxval(abs(dense(A) - Amat)) < n**2 * epsilon(1.0_dp), &
                 "SymTridiagonal eigh failed.")
      return
   end subroutine test_eigh
end module
