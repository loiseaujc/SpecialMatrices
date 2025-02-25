module test_poisson2D
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

   integer, parameter :: nx = 16, ny = 16, n = nx*ny

   public :: collect_poisson2D_testsuite

contains

   !-----------------------------------
   !-----     poisson2D MATRICES     -----
   !-----------------------------------
   subroutine collect_poisson2D_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Poisson2D trace", test_trace), &
                  ! new_unittest("Poisson2D determinant", test_det), &
                  new_unittest("Poisson2D matmul", test_matmul), &
                  ! new_unittest("Poisson2D linear solver", test_solve), &
                  new_unittest("Poisson2D eigenvalue decomposition", test_eigh) &
                  ]
      return
   end subroutine collect_poisson2D_testsuite

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Poisson2D) :: A

      ! Initialize matrix.
      A = Poisson2D(nx, ny)
      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Poisson2D matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Poisson2D matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(Poisson2D) :: A

      ! Initialize matrix.
      A = Poisson2D(nx, ny)

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
                    "Poisson2D solve with a single rhs failed.")
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
                    "Poisson2D solve with multiple rhs failed.")
      end block

      return
   end subroutine test_solve

   subroutine test_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(Poisson2D) :: A

      ! Initialize matrix.
      A = Poisson2D(nx, ny)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "Poisson2D trace failed.")
      return
   end subroutine test_trace

   subroutine test_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(Poisson2D) :: A

      ! Initialize matrix.
      A = Poisson2D(nx, ny)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "Poisson2D det failed.")
      return
   end subroutine test_det

   subroutine test_eigh(error)
      type(error_type), allocatable, intent(out) :: error
      type(Poisson2D) :: A
      real(dp), allocatable :: lambda(:), vectors(:, :), Amat(:, :)
      real(dp), allocatable :: lambda_stdlib(:)
      integer :: i, j

      ! Initialize matrix.
      A = Poisson2D(nx, ny)

      ! Compute eigenvalues.
      lambda = eigvalsh(A); lambda_stdlib = eigvalsh(dense(A))
      ! Check error.
      call check(error, all_close(lambda, lambda_stdlib), &
                 "Poisson2D eigvalsh failed.")
      if (allocated(error)) return

      ! Compute eigenvalues and eigenvectors.
      call eigh(A, lambda, vectors)
      ! Check error.
      Amat = matmul(vectors, matmul(diag(lambda), transpose(vectors)))
      call check(error, maxval(abs(dense(A) - Amat)) < n**2*epsilon(1.0_dp), &
                 "Poisson2D eigh failed.")
      return
   end subroutine test_eigh
end module
