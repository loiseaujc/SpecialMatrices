module test_hankel
   ! Fortran standard library.
   use stdlib_math, only: is_close, all_close
   use stdlib_sorting, only: sort_index
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg, only: norm, diag, svdvals, eye, mnorm
   ! Testdrive.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   ! SpecialMatrices
   use SpecialMatrices, only: hankel, dense, transpose, matmul, operator(*), &
                              svd, svdvals 
   implicit none(type, external)
   private

   public :: collect_hankel_testsuite
contains

   !-------------------------------------
   !-----      HANKEL MATRICES      -----
   !-------------------------------------

   subroutine collect_hankel_testsuite(testsuite)
      implicit none(type, external)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Hankel matrix transpose", test_transpose), &
                  new_unittest("Hankel scalar multiplication", test_scalar_multiplication), &
                  new_unittest("Hankel matmul", test_matmul), &
                  new_unittest("Hankel svdvals", test_svdvals), &
                  new_unittest("Hankel svd", test_svd) &
                  ]
      return
   end subroutine collect_hankel_testsuite

   subroutine test_transpose(error)
      implicit none(type, external)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: m = 2, n = 3
      type(hankel) :: H
      real(dp) :: v(m + n - 1), A(m, n)

      ! Initialize matrix.
      A(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]
      A(2, :) = [2.0_dp, 3.0_dp, 4.0_dp]

      ! Initialize vector.
      v = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]

      ! Initialize Hankel matrix.
      H = Hankel(v, m, n)

      ! Check error.
      call check(error, all_close(dense(transpose(H)), transpose(A)), &
                 "Transposition failed.")
   end subroutine test_transpose

   subroutine test_scalar_multiplication(error)
      implicit none(type, external)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: m = 128, n = 256
      type(hankel) :: A, B
      real(dp), allocatable :: v(:)
      real(dp) :: alpha

      ! Initialize matrix.
      allocate (v(m + n - 1)); call random_number(v)
      A = Hankel(v, m, n); call random_number(alpha)

      ! Scalar-matrix multiplication.
      B = alpha*A

      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "alpha*Hankel failed.")
      if (allocated(error)) return

      ! Matrix-scalar multipliation.
      B = A*alpha
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "Hankel*alpha failed.")

      return
   end subroutine test_scalar_multiplication

   subroutine test_matmul(error)
      implicit none(type, external)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: m = 8, n = 6
      type(hankel) :: A
      real(dp), allocatable :: v(:)

      ! Initialize matrix.
      allocate (v(m + n - 1)); call random_number(v)
      A = Hankel(v, m, n)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n)); call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Hankel matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Hankel matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_svdvals(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: m = 8, n = 6
      type(hankel) :: A
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: s(:), s_stdlib(:)

      ! Initialize Matrix.
      allocate (v(m + n -1)); call random_number(v)
      A = Hankel(v, m, n)
      ! Compute singular values.
      s = svdvals(A); s_stdlib = svdvals(dense(A))
      ! Check error.
      call check(error, all_close(s, s_stdlib), &
                 "hankel svdvals failed.")
      return
   end subroutine test_svdvals

   subroutine test_svd(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: m = 8, n = 6
      integer :: i
      type(hankel) :: A
      real(dp), allocatable :: v(:), Amat(:, :)
      real(dp), allocatable :: u(:, :), s(:), vt(:, :)
      real(dp), allocatable :: sigma(:, :)

      !Initialize Matrix.
      allocate(v(m + n -1)); call random_number(v)
      A = Hankel(v, m, n)

      ! Compute singular value decomosition.
      allocate(s(min(m, n)), u(m, m), vt(n, n))
      call svd(A, s, u, vt)

      ! Check orthogonality of the left singular vectors.
      block
         real(dp), allocatable :: G(:, :)
         G = matmul(transpose(u), u)
         call check(error, norm(G - eye(m, mold=1.0_dp), "inf") < 1e-10_dp, &
                    "Orthogonality of the left singular vectors failed.")
      end block

      ! Check orthogonality of the right singular vectors.
      block
         real(dp), allocatable :: G(:, :)
         G = matmul(vt, transpose(vt))
         call check(error, norm(G - eye(n, mold=1.0_dp), "inf") < 1e-10_dp, &
                    "Orthogonality of the right singular vectors failed.")
      end block
       
      ! Build Sigma
      allocate(sigma(m, n)); sigma = 0.0_dp
      do i = 1, min(m, n)
         sigma(i, i) = s(i)
      end do

      ! Check error.
      allocate (Amat(m, n)); Amat = 0.0_dp
      Amat = matmul(u, matmul(sigma, vt)) ! s = mxn or full_matrices=.false.
      call check(error, mnorm(dense(A) - Amat, 2) < 1e-8_dp, &
                 "hankel svd failed.")
      return
   end subroutine test_svd


end module test_hankel
