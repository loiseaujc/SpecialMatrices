module test_bidiagonal
   ! Fortran standard library.
   use stdlib_math, only: is_close, all_close
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg, only: diag, det, trace, inv, solve, svdvals, eigvals, eig
   ! Testdrive.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   ! SpecialMatrices
   use SpecialMatrices
   implicit none
   private

   integer, parameter :: n = 512

   public :: collect_bidiagonal_testsuite
contains

   !---------------------------------------
   !-----     BIDIAGONAL MATRICES     -----
   !---------------------------------------
   subroutine collect_bidiagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Bidiagonal scalar multiplication", test_scalar_multiplication), &
                  new_unittest("Bidiagonal trace", test_trace), &
                  new_unittest("Bidiagonal determinant", test_det), &
                  new_unittest("Bidiagonal matmul", test_matmul), &
                  new_unittest("Bidiagonal linear solver", test_solve), &
                  ! new_unittest("Bidiagonal eigenvalue decomposition", test_eig), &
                  new_unittest("Bidiagonal singular value decomposition", test_svd) &
]
      return
   end subroutine collect_bidiagonal_testsuite

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1))
      call random_number(dv); call random_number(ev)
      A = Bidiagonal(dv, ev, which="L")

      ! Matrix-vector product with A lower bidiagonal.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Lower-bidiagonal matrix-vector product failed.")
         if (allocated(error)) return
         ! Matrix-vector product wih A upper bidiaognal.
         A = Bidiagonal(dv, ev, which="U")
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Upper-bidiagonal matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product with A upper bidiagonal.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Upper-bidiagonal matrix-matrix product failed.")
         if (allocated(error)) return
         ! Matrix-matrix product with A lower bidiagonal.
         A = Bidiagonal(dv, ev, which="L")
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Lower-bidiagonal matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A
      real(dp), allocatable :: ev(:), dv(:)

      ! Initialize matrix.
      allocate (ev(n - 1), dv(n))
      call random_number(ev); call random_number(dv)
      A = Bidiagonal(dv, ev, which="L")

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
                    "Lower-bidiagonal solve with a single rhs failed.")
         if (allocated(error)) return

         A = Bidiagonal(dv, ev, which="U")
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
                    "Upper-bidiagonal solve with multiple rhs failed.")
         if (allocated(error)) return

         A = Bidiagonal(dv, ev, which="L")
         ! Solve with SpecialMatrices.
         x = solve(A, b)
         ! Solve with stdlib (dense).
         x_stdlib = solve(dense(A), b)
         ! Check error.
         call check(error, all_close(x, x_stdlib), &
                    "Lower-bidiagonal solve with multiple rhs failed.")
      end block

      return
   end subroutine test_solve

   subroutine test_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n-1)); call random_number(dv); call random_number(ev)
      A = Bidiagonal(dv, ev)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
      "Bidiagonal trace failed.")
      
   end subroutine test_trace

   subroutine test_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n-1)); call random_number(dv); call random_number(ev)
      A = Bidiagonal(dv, ev)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
      "Bidiagonal det failed.")
      
   end subroutine test_det

   subroutine test_scalar_multiplication(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A, B
      real(dp), allocatable :: dv(:), ev(:)
      real(dp) :: alpha

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1)); call random_number(dv); call random_number(ev)
      A = Bidiagonal(dv, ev)
      call random_number(alpha)

      ! Scalar-matrix multiplication.
      B = alpha*A
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "alpha*Bidiagonal failed.")
      if (allocated(error)) return

      ! Matrix-scalar multipliation.
      B = A*alpha
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "Bidiagonal*alpha failed.")
      return
   end subroutine test_scalar_multiplication

   subroutine test_eig(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A
      real(dp), allocatable :: ev(:), dv(:)
      complex(dp), allocatable :: lambda(:), right(:, :), left(:, :)
      complex(dp), allocatable :: Amat(:, :), diag_a(:)
      integer :: i

      ! Initialize matrix.
      allocate (ev(n - 1), dv(n))
      call random_number(dv); call random_number(ev)
      A = Bidiagonal(dv, ev, which="U")

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
                 "Bidiagonal eig failed.")
      return
   end subroutine test_eig

   subroutine test_svd(error)
      type(error_type), allocatable, intent(out) :: error
      type(Bidiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)
      real(dp), allocatable :: s(:), u(:, :), vt(:, :)
      real(dp), allocatable :: Amat(:, :)

      ! Initialize matrix.
      allocate (dv(n), ev(n - 1))
      call random_number(dv); call random_number(ev)
      A = Bidiagonal(dv, ev)

      ! Compute eigendecomposition.
      allocate (s(n), u(n, n), vt(n, n))
      call svd(A, s, u, vt)

     ! Check error.
      Amat = matmul(u, matmul(diag(s), vt))
      call check(error, maxval(abs(dense(A) - Amat)) < 10*n**2*epsilon(1.0_dp), &
                 "Bidiagonal svd failed.")

      return
   end subroutine test_svd
end module
