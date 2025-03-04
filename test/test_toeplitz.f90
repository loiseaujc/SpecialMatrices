module test_toeplitz
   ! Fortran standard library.
   use stdlib_math, only: is_close, all_close
   use stdlib_sorting, only: sort_index
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg, only: norm
   ! Testdrive.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   ! SpecialMatrices
   use SpecialMatrices
   implicit none
   private

   integer, parameter :: m = 128, n = 128

   public :: collect_toeplitz_testsuite
contains

   !-------------------------------------
   !-----     TOEPLITZ MATRICES     -----
   !-------------------------------------

   subroutine collect_toeplitz_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Toeplitz scalar multiplication", test_scalar_multiplication), &
                  new_unittest("Toeplitz matmul", test_matmul), &
                  new_unittest("Toeplitz linear solver", test_solve) &
                  ]
      return
   end subroutine collect_toeplitz_testsuite

   subroutine test_scalar_multiplication(error)
      type(error_type), allocatable, intent(out) :: error
      type(toeplitz) :: A, B
      real(dp), allocatable :: vc(:), vr(:)
      real(dp) :: alpha

      ! Initialize matrix.
      allocate(vc(m)) ; call random_number(vc)
      allocate(vr(n)) ; call random_number(vr)
      A = Toeplitz(vc, vr) ; call random_number(alpha)

      ! Scalar-matrix multiplication.
      B = alpha*A
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "alpha*toeplitz failed.")
      if (allocated(error)) return

      ! Matrix-scalar multipliation.
      B = A*alpha
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "toeplitz*alpha failed.")
      return
   end subroutine test_scalar_multiplication

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Toeplitz) :: A
      real(dp), allocatable :: vc(:), vr(:)

      ! Initialize matrix.
      allocate(vc(m)) ; call random_number(vc)
      allocate(vr(n)) ; call random_number(vr)
      A = Toeplitz(vc, vr)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n)); call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Toeplitz matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Toeplitz matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(toeplitz) :: A
      real(dp), allocatable :: vr(:), vc(:)
      integer(ilp) :: i

      ! Initialize matrix.
      vr = [(1.0_dp / (i+1), i=1, n)]
      vc = [(1.0_dp / (i+1), i=1, n)]
      A = Toeplitz(vc, vr)

      ! Solve with a single right-hand side vector.
      block
         real(dp), allocatable :: x(:), b(:)
         allocate (b(n))
         ! Random rhs.
         call random_number(b) ; b = b / norm(b, 2)
         ! Solve with SpecialMatrices.
         x = solve(A, b)
         ! Check error.
         call check(error, norm(matmul(A, x) - b, 2) <= 1e-8_dp, &
                    "toeplitz solve with a single rhs failed.")
         if (allocated(error)) return
      end block

      ! Solve with multiple right-hand side vectors.
      block
         real(dp), allocatable :: x(:, :), b(:, :)
         integer(ilp) :: i
         allocate (b(n, n))
         ! Random rhs.
         call random_number(b)
         do i = 1, n
            b(:, i) = b(:, i) / norm(b, 2)
         enddo
         ! Solve with SpecialMatrices.
         x = solve(A, b)
         ! Check error.
         do i = 1, n
            call check(error, norm(matmul(A, x(:, i)) - b(:, i), 2) <= 1e-8_dp, &
                      "toeplitz solve with multiple rhs failed.")
            if (allocated(error)) return
         enddo
      end block

      return
   end subroutine test_solve

end module
