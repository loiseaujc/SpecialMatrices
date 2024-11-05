module TestTridiag
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
   public :: collect_bidiagonal_testsuite
   public :: collect_tridiagonal_testsuite
   public :: collect_symtridiagonal_testsuite

contains

   !-------------------------------------
   !-----     DIAGONAL MATRICES     -----
   !-------------------------------------

   subroutine collect_diagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Diagonal scalar multiplication", test_diagonal_scalar_multiplication), &
                  new_unittest("Diagonal trace", test_diagonal_trace), &
                  new_unittest("Diagonal determinant", test_diagonal_det), &
                  new_unittest("Diagonal inverse", test_diagonal_inv), &
                  new_unittest("Diagonal matmul", test_diagonal_matmul), &
                  new_unittest("Diagonal linear solver", test_diagonal_solve), &
                  new_unittest("Diagonal svdvals", test_diagonal_svdvals), &
                  new_unittest("Diagonal svd", test_diagonal_svd), &
                  new_unittest("Diagonal eigvalsh", test_diagonal_eigvalsh), &
                  new_unittest("Diagonal eigh", test_diagonal_eigh) &
                  ]
      return
   end subroutine collect_diagonal_testsuite

   subroutine test_diagonal_scalar_multiplication(error)
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
   end subroutine test_diagonal_scalar_multiplication

   subroutine test_diagonal_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "Diagonal trace failed.")
      return
   end subroutine test_diagonal_trace

   subroutine test_diagonal_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "Diagonal det failed.")
      return
   end subroutine test_diagonal_det

   subroutine test_diagonal_inv(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Intialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Compare against stdlib_linalg implementation.
      call check(error, all_close(inv(dense(A)), dense(inv(A))), &
                 "Diagonal inv failed.")
      return
   end subroutine test_diagonal_inv

   subroutine test_diagonal_matmul(error)
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
   end subroutine test_diagonal_matmul

   subroutine test_diagonal_solve(error)
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
   end subroutine test_diagonal_solve

  subroutine test_diagonal_svdvals(error)
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
   end subroutine test_diagonal_svdvals

   subroutine test_diagonal_svd(error)
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
   end subroutine test_diagonal_svd

   subroutine test_diagonal_eigvalsh(error)
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
   end subroutine test_diagonal_eigvalsh

   subroutine test_diagonal_eigh(error)
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
   end subroutine test_diagonal_eigh

   !---------------------------------------
   !-----     BIDIAGONAL MATRICES     -----
   !---------------------------------------
   subroutine collect_bidiagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Bidiagonal matmul", test_bidiagonal_matmul), &
                  new_unittest("Bidiagonal linear solver", test_bidiagonal_solve) &
                  ]
      return
   end subroutine collect_bidiagonal_testsuite

   subroutine test_bidiagonal_matmul(error)
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
   end subroutine test_Bidiagonal_matmul

   subroutine test_bidiagonal_solve(error)
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
   end subroutine test_bidiagonal_solve

   !----------------------------------------
   !-----     TRIDIAGONAL MATRICES     -----
   !----------------------------------------

   subroutine collect_tridiagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("Tridiagonal matmul", test_tridiagonal_matmul), &
                  new_unittest("Tridiagonal linear solver", test_tridiagonal_solve) &
                  ]
      return
   end subroutine collect_tridiagonal_testsuite

   subroutine test_tridiagonal_matmul(error)
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
   end subroutine test_tridiagonal_matmul

   subroutine test_tridiagonal_solve(error)
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
         x = solve(A, b)
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
   subroutine collect_symtridiagonal_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("SymTridiagonal scalar multiplication", test_symtridiagonal_scalar_multiplication), &
                  new_unittest("SymTridiagonal trace", test_symtridiagonal_trace), &
                  new_unittest("SymTridiagonal determinant", test_symtridiagonal_det), &
                  new_unittest("SymTridiagonal matmul", test_symtridiagonal_matmul), &
                  new_unittest("SymTridiagonal linear solver", test_symtridiagonal_solve), &
                  new_unittest("Pos. def. SymTridiagonal linear solver", test_posdef_symtridiagonal_solve) &
                  ]
      return
   end subroutine collect_symtridiagonal_testsuite

   subroutine test_symtridiagonal_matmul(error)
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
   end subroutine test_symtridiagonal_matmul

   subroutine test_symtridiagonal_solve(error)
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
   end subroutine test_symtridiagonal_solve

   subroutine test_posdef_symtridiagonal_solve(error)
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
   end subroutine test_posdef_symtridiagonal_solve

   subroutine test_symtridiagonal_scalar_multiplication(error)
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
   end subroutine test_symtridiagonal_scalar_multiplication

   subroutine test_symtridiagonal_trace(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n-1)); call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(trace(A), trace(dense(A))), &
                 "SymTridiagonal trace failed.")
      return
   end subroutine test_symtridiagonal_trace

   subroutine test_symtridiagonal_det(error)
      type(error_type), allocatable, intent(out) :: error
      type(SymTridiagonal) :: A
      real(dp), allocatable :: dv(:), ev(:)

      ! Initialize matrix.
      allocate (dv(n), ev(n-1)); call random_number(dv); call random_number(ev)
      A = SymTridiagonal(dv, ev)

      ! Compare against stdlib_linalg implementation.
      call check(error, is_close(det(A), det(dense(A))), &
                 "SymTridiagonal det failed.")
      return
   end subroutine test_symtridiagonal_det
end module
