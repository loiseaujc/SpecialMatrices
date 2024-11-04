module TestTridiag
   ! Fortran standard library.
   use stdlib_math, only: is_close, all_close
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg, only: diag, det, trace, inv, solve, svdvals
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
                  new_unittest("Diagonal trace", test_diagonal_trace), &
                  new_unittest("Diagonal determinant", test_diagonal_det), &
                  new_unittest("Diagonal inverse", test_diagonal_inv), &
                  new_unittest("Diagonal matmul", test_diagonal_matmul), &
                  new_unittest("Diagonal in-place matmul", test_diagonal_spmv_ip), &
                  new_unittest("Diagonal linear solver", test_diagonal_solve), &
                  new_unittest("Diagonal in-place linear solver", test_diagonal_solve_ip), &
                  new_unittest("Diagonal svdvals", test_diagonal_svdvals), &
                  new_unittest("Diagonal svd", test_diagonal_svd) &
                  ]
      return
   end subroutine collect_diagonal_testsuite

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
      call check(error, all_close(inv(dense(A)), inv(A)), &
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

   subroutine test_diagonal_spmv_ip(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n), y(n)); call random_number(x)
         call spmv_ip(y, A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "in-place Diagonal matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n), y(n, n))
         call random_number(x)
         call spmv_ip(y, A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "in-place Diagonal matrix-matrix product failed.")
      end block
      return
   end subroutine test_diagonal_spmv_ip

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

   subroutine test_diagonal_solve_ip(error)
      type(error_type), allocatable, intent(out) :: error
      type(Diagonal) :: A
      real(dp), allocatable :: dv(:)

      ! Initialize matrix.
      allocate (dv(n)); call random_number(dv); A = Diagonal(dv)

      ! Solve with a single right-hand side vector.
      block
         real(dp), allocatable :: x(:), x_stdlib(:), b(:)
         allocate (b(n), x(n))
         ! Random rhs.
         call random_number(b)
         ! Solve with SpecialMatrices.
         call solve_ip(x, A, b)
         ! Solve with stdlib (dense).
         x_stdlib = solve(dense(A), b)
         ! Check error.
         call check(error, all_close(x, x_stdlib), &
                    "in-place Diagonal solve with a single rhs failed.")
         if (allocated(error)) return
      end block

      ! Solve with multiple right-hand side vectors.
      block
         real(dp), allocatable :: x(:, :), x_stdlib(:, :), b(:, :)
         allocate (b(n, n), x(n, n))
         ! Random rhs.
         call random_number(b)
         ! Solve with SpecialMatrices.
         call solve_ip(x, A, b)
         ! Solve with stdlib (dense).
         x_stdlib = solve(dense(A), b)
         ! Check error.
         call check(error, all_close(x, x_stdlib), &
                    "in-place Diagonal solve with multiple rhs failed.")
      end block

      return
   end subroutine test_diagonal_solve_ip

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
                  new_unittest("SymTridiagonal matmul", test_symtridiagonal_matmul), &
                  new_unittest("SymTridiagonal linear solver", test_symtridiagonal_solve) &
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
         allocate (x(n, n), y(n, n), y_dense(n, n))
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

end module
