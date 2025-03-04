module test_diagonal
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

   public :: collect_diagonal_testsuite
contains

   !-------------------------------------
   !-----     DIAGONAL MATRICES     -----
   !-------------------------------------

   subroutine collect_diagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Diagonal scalar multiplication", test_scalar_multiplication), &
                  new_unittest("Diagonal trace", test_trace), &
                  new_unittest("Diagonal determinant", test_det), &
                  new_unittest("Diagonal inverse", test_inv), &
                  new_unittest("Diagonal matmul", test_matmul), &
                  new_unittest("Diagonal linear solver", test_solve), &
                  new_unittest("Diagonal svdvals", test_svdvals), &
                  new_unittest("Diagonal svd", test_svd), &
                  new_unittest("Diagonal eigvalsh", test_eigvalsh), &
                  new_unittest("Diagonal eigh", test_eigh) &
                  ]
      return
   end subroutine collect_diagonal_testsuite

   subroutine test_scalar_multiplication(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A, B
      real(dp), allocatable :: dv(:)
      real(dp) :: alpha

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)
      call random_number(alpha)

      ! Scalar-matrix multiplication.
      B = alpha*A
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "alpha*Diagonal failed.")
      if (allocated(error)) return

      ! Matrix-scalar multipliation.
      B = A*alpha
      ! Check error.
      call check(error, all_close(alpha*dense(A), dense(B)), &
                 "Diagonal*alpha failed.")
      return
   end subroutine test_scalar_multiplication

   subroutine test_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "Diagonal trace failed.")
      return
   end subroutine test_trace

   subroutine test_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "Diagonal det failed.")
      return
   end subroutine test_det

   subroutine test_inv(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Intialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Compare against stdlib_linalg implementation.
      call check(error, all_close(inv(dense(A)), dense(inv(A))), &
                 "Diagonal inv failed.")
      return
   end subroutine test_inv

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n)); call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Diagonal matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Diagonal matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   subroutine test_solve(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Solve with a single right-hand side vector.
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
                    "Diagonal solve with a single rhs failed.")
         if (allocated(error)) return
      end block

      ! Solve with multiple right-hand side vectors.
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
                    "Diagonal solve with multiple rhs failed.")
      end block

      return
   end subroutine test_solve

   subroutine test_svdvals(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)
      real(dp), allocatable :: s(:), s_stdlib(:)

      ! Initialize array.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)
      ! Compute singular values.
      s = svdvals(A); s_stdlib = svdvals(dense(A))
      ! Check error.
      call check(error, all_close(s, s_stdlib), &
                 "Diagonal svdvals failed.")
      return
   end subroutine test_svdvals

   subroutine test_svd(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:), Amat(:, :)
      real(dp), allocatable :: u(:, :), s(:), vt(:, :)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)
      ! Compute singular value decomposition.
      call svd(A, u, s, vt)
      ! Check error.
      allocate (Amat(n, n)); Amat = 0.0_dp
      Amat = matmul(u, matmul(diag(s), vt))
      call check(error, all_close(dense(A), Amat), &
                 "Diagonal svd failed.")
      return
   end subroutine test_svd

   subroutine test_eigvalsh(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)
      real(dp), allocatable :: lambda(:), lambda_stdlib(:)

      ! Initialize array.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)
      ! Compute singular values.
      lambda = eigvalsh(A); lambda_stdlib = eigvalsh(dense(A))
      ! Check error.
      call check(error, all_close(lambda, lambda_stdlib), &
                 "Diagonal eigvalsh failed.")
      return
   end subroutine test_eigvalsh

   subroutine test_eigh(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:), Amat(:, :)
      real(dp), allocatable :: lambda(:), vectors(:, :)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)
      ! Compute singular value decomposition.
      call eigh(A, lambda, vectors)
      ! Check error.
      allocate (Amat(n, n)); Amat = 0.0_dp
      Amat = matmul(vectors, matmul(diag(lambda), transpose(vectors)))
      call check(error, all_close(dense(A), Amat), &
                 "Diagonal eigh failed.")
      return
   end subroutine test_eigh

end module
