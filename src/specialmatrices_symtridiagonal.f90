submodule(SpecialMatrices_Tridiagonal) SymmetricTridiagonal
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_symtridiag(n) result(A)
      ! Dimension of the matrix.
      integer(ilp), intent(in) :: n
      ! Output matrix.
      type(SymTridiagonal) :: A
      A%n = n; allocate(A%dv(n), A%ev(n-1))
      A%dv = 0.0_dp; A%ev = 0.0_dp
      return
   end function initialize_symtridiag

   pure module function construct_symtridiag(dv, ev) result(A)
      ! Diagonals elements.
      real(dp), intent(in) :: dv(:), ev(:)
      ! Output matrix.
      type(SymTridiagonal) :: A
      A%n = size(dv); A%dv = dv; A%ev = ev
      return
   end function construct_symtridiag

   pure module function construct_constant_symtridiag(d, e, n) result(A)
      ! Diagonal elements.
      real(dp), intent(in) :: d, e
      ! Dimension of the matrix.
      integer(ilp), intent(in) :: n
      ! Output matrix.
      type(SymTridiagonal) :: A
      integer i
      A%n = n; A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n-1)]
      return
   end function construct_constant_symtridiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function symtridiag_spmv(A, x) result(y)
      ! Input matrix.
      type(SymTridiagonal), intent(in) :: A
      ! Input vector.
      real(dp), intent(in) :: x(:)
      ! Output vector.
      real(dp) :: y(size(x))
      integer :: i, n
      n = size(x); y(1) = A%dv(1)*x(1) + A%ev(1)*x(2)
      do concurrent (i=2:n-1)
         y(i) = A%ev(i-1)*x(i-1) + A%dv(i)*x(i) + A%ev(i)*x(i+1)
      enddo
      y(n) = A%dv(n)*x(n) + A%ev(n-1)*x(n-1)
      return
   end function symtridiag_spmv

   pure module function symtridiag_multi_spmv(A, X) result(Y)
      ! Coefficient matrix.
      type(SymTridiagonal), intent(in) :: A
      ! Input vectors.
      real(dp), intent(in) :: X(:, :)
      ! Output vectors.
      real(dp) :: Y(size(X, 1), size(X, 2))
      integer(ilp) :: i
      do concurrent(i=1:size(X,2))
         Y(:, i) = symtridiag_spmv(A, X(:, i))
      enddo
   end function symtridiag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   pure module function symtridiag_solve(A, b) result(x)
      ! Coefficient matrix.
      type(SymTridiagonal), intent(in) :: A
      ! Right-hand side vector.
      real(dp), intent(in) :: b(:)
      ! Solution vector.
      real(dp) :: x(size(b))
      integer :: i, n, nrhs, info
      real(dp) :: dl(A%n - 1), d(A%n), du(A%n - 1), b_(A%n, 1)

      ! Initialize arrays.
      n = A%n; dl = A%ev; d = A%dv; du = A%ev; b_(:, 1) = b ; nrhs = 1

      ! Solve the system.
      call gtsv(n, nrhs, dl, d, du, b_, n, info)

      ! Return the solution.
      x = b_(:, 1)

      return
   end function symtridiag_solve

   pure module function symtridiag_multi_solve(A, B) result(X)
      ! Coefficient matrix.
      type(SymTridiagonal), intent(in) :: A
      ! Right-hand side vectors.
      real(dp), intent(in) :: B(:, :)
      ! Solution vectors.
      real(dp) :: X(size(B, 1), size(B, 2))
      integer :: i, n, nrhs, info
      real(dp) :: dl(A%n - 1), d(A%n), du(A%n - 1)

      ! Initialize arrays.
      n = A%n; dl = A%ev; d = A%dv; du = A%ev; nrhs = size(B, 2); X = B

      ! Solve the systems.
      call gtsv(n, nrhs, dl, d, du, X, n, info)

   end function symtridiag_multi_solve

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   pure module function symtridiag_to_dense(A) result(B)
      ! Input tridiagonal matrix.
      type(SymTridiagonal), intent(in) :: A
      ! Output dense matrix.
      real(dp) :: B(A%n, A%n)
      integer :: i, n
      n = A%n; B = 0.0_dp
      B(1, 1) = A%dv(1); B(1, 2) = A%ev(1)
      do concurrent(i=2:n - 1)
         B(i, i - 1) = A%ev(i - 1)
         B(i, i) = A%dv(i)
         B(i, i + 1) = A%ev(i)
      end do
      B(n, n - 1) = A%ev(n - 1); B(n, n) = A%dv(n)
      return
   end function symtridiag_to_dense

   pure module function symtridiag_transpose(A) result(B)
      ! Input tridiagonal matrix.
      type(SymTridiagonal), intent(in) :: A
      ! Transposed matrix.
      type(SymTridiagonal) :: B
      B = A
   end function symtridiag_transpose

end submodule SymmetricTridiagonal
