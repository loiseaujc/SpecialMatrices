module test_strang
   ! Fortran standard library.
   use stdlib_io_npy, only: save_npy
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

   public :: collect_strang_testsuite

contains

   !-----------------------------------
   !-----     STRANG MATRICES     -----
   !-----------------------------------
   subroutine collect_strang_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Strang trace", test_trace), &
                  new_unittest("Strang determinant", test_det), &
                  new_unittest("Strang matmul", test_matmul), &
                  new_unittest("Strang linear solver", test_solve), &
                  new_unittest("Strang eigenvalue decomposition", test_eigh) &
                  ]
      return
   end subroutine collect_strang_testsuite

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Strang) :: A

      ! Initialize matrix.
      A = Strang(n)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Strang matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Strang matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(Strang) :: A

      ! Initialize matrix.
      A = Strang(n)

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
                    "Strang solve with a single rhs failed.")
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
                    "Strang solve with multiple rhs failed.")
      end block

      return
   end subroutine test_solve

   subroutine test_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(Strang) :: A

      ! Initialize matrix.
      A = Strang(n)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "Strang trace failed.")
      return
   end subroutine test_trace

   subroutine test_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(Strang) :: A

      ! Initialize matrix.
      A = Strang(n)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "Strang det failed.")
      return
   end subroutine test_det

   subroutine test_eigh(error)
      type(error_type), allocatable, intent(out) :: error
      type(Strang) :: A
      real(dp), allocatable :: lambda(:), vectors(:, :), Amat(:, :)
      real(dp), allocatable :: lambda_stdlib(:)
      integer :: i, j

      ! Initialize matrix.
      A = Strang(n)

      ! Compute eigenvalues.
      lambda = eigvalsh(A); lambda_stdlib = eigvalsh(dense(A))

      ! Check error.
      call check(error, all_close(lambda, lambda_stdlib), &
                 "Strang eigvalsh failed.")
      if (allocated(error)) return

      ! Compute eigenvalues and eigenvectors.
      call eigh(A, lambda, vectors)
      ! Check error.
      Amat = matmul(vectors, matmul(diag(lambda), transpose(vectors)))
      call check(error, maxval(abs(dense(A) - Amat)) < n**2*epsilon(1.0_dp), &
                 "Strang eigh failed.")
      return
   end subroutine test_eigh
end module
