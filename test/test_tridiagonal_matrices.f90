module test_tridiagonal
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

   public :: collect_tridiagonal_testsuite

contains

   !----------------------------------------
   !-----     TRIDIAGONAL MATRICES     -----
   !----------------------------------------

   subroutine collect_tridiagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Tridiagonal scalar multiplication", test_scalar_multiplication), &
                  new_unittest("Tridiagonal trace", test_trace), &
                  new_unittest("Tridiagonal determinant", test_det), &
                  new_unittest("Tridiagonal matmul", test_matmul), &
                  new_unittest("Tridiagonal linear solver", test_solve), &
                  new_unittest("Tridiagonal eigenvalue decomposition", test_eig), &
                  new_unittest("Tridiagonal singular value decomposition", test_svd) &
                  ]
      return
   end subroutine collect_tridiagonal_testsuite

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A
      real(dp), allocatable :: dl(:), dv(:), du(:)

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = Tridiagonal(dl, dv, du)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Tridiagonal matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Tridiagonal matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A
      real(dp), allocatable :: dl(:), dv(:), du(:)

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = Tridiagonal(dl, dv, du)

      ! Solve with a singe right-hand side vector.
      block
         real(dp), allocatable :: x(:), x_stdlib(:), b(:)
         allocate (b(n))
         ! Random rhs.
         call random_number(b)
         ! Solve with SpecialMatrices.
         x = solve(A, b, refine=.true.)
         ! Solve with stdlib (dense).
         x_stdlib = solve(dense(A), b)
         ! Check error.
         call check(error, all_close(x, x_stdlib), &
                    "Tridiagonal solve with a single rhs failed.")
         if (allocated(error)) return
      end block

      block
         real(dp), allocatable :: x(:, :), x_stdlib(:, :), b(:, :)
         allocate (b(n, n))
         ! Random rhs.
         call random_number(b)
         ! Solve with SpecialMatrices.
         x = solve(A, b, refine=.true.)
         ! Solve with stdlib (dense).
         x_stdlib = solve(dense(A), b)
         ! Check error.
         call check(error, all_close(x, x_stdlib), &
                    "Tridiagonal solve with multiple rhs failed.")
      end block

      return
   end subroutine test_solve

   subroutine test_scalar_multiplication(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A, B
      real(dp), allocatable :: dl(:), dv(:), du(:)
      real(dp) :: alpha

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = TriDiagonal(dl, dv, du)
      call random_number(alpha)

      ! Scalar-matrix multiplication.
      B = alpha*A
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "alpha*Tridiagonal failed.")
      if (allocated(error)) return

      ! Matrix-scalar multipliation.
      B = A*alpha
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "Tridiagonal*alpha failed.")
      return
   end subroutine test_scalar_multiplication

   subroutine test_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A
      real(dp), allocatable :: dl(:), dv(:), du(:)

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = Tridiagonal(dl, dv, du)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "Tridiagonal trace failed.")
      return
   end subroutine test_trace

   subroutine test_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A
      real(dp), allocatable :: dl(:), dv(:), du(:)

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = Tridiagonal(dl, dv, du)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "SymTridiagonal det failed.")
      return
   end subroutine test_det

   subroutine test_eig(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A
      real(dp), allocatable :: dl(:), dv(:), du(:)
      complex(dp), allocatable :: lambda(:), right(:, :), left(:, :)
      complex(dp), allocatable :: Amat(:, :), diag_a(:)
      integer :: i

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = Tridiagonal(dl, dv, du)

      ! Compute eigendecomposition.
      allocate (lambda(n), left(n, n), right(n, n))
      call eig(A, lambda, left=left, right=right)

      ! Normalize eigenvectors.
      do i = 1, n
         right(:, i) = right(:, i)/dot_product(right(:, i), left(:, i))
      end do

      ! Check error.
      Amat = matmul(right, matmul(diag(lambda), conjg(transpose(left))))
      call check(error, maxval(abs(dense(A) - Amat%re)) < 10*n**2*epsilon(1.0_dp), &
                 "Tridiagonal eig failed.")

      return
   end subroutine test_eig

   subroutine test_svd(error)
      type(error_type), allocatable, intent(out) :: error
      type(Tridiagonal) :: A
      real(dp), allocatable :: dl(:), dv(:), du(:)
      real(dp), allocatable :: s(:), u(:, :), vt(:, :)
      real(dp), allocatable :: Amat(:, :)

      ! Initialize matrix.
      allocate (dl(n - 1), dv(n), du(n - 1))
      call random_number(dl); call random_number(dv); call random_number(du)
      A = Tridiagonal(dl, dv, du)

      ! Compute eigendecomposition.
      allocate (s(n), u(n, n), vt(n, n))
      call svd(A, s, u, vt)

     ! Check error.
      Amat = matmul(u, matmul(diag(s), vt))
      call check(error, maxval(abs(dense(A) - Amat)) < 10*n**2*epsilon(1.0_dp), &
                 "Tridiagonal svd failed.")

      return
   end subroutine test_svd
end module
