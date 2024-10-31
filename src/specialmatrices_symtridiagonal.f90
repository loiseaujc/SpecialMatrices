submodule(SpecialMatrices_Tridiagonal) SymmetricTridiagonal
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_symtridiag(n) result(A)
      !! Utility function to create a `SymTridiagonal` matrix filled with zeros.
      integer(ilp), intent(in) :: n
      !! Dimension of the matrix.
      type(SymTridiagonal) :: A
      !! Corresponding symmetric tridiagonal matrix.
      A%n = n; allocate(A%dv(n), A%ev(n-1))
      A%dv = 0.0_dp; A%ev = 0.0_dp
      return
   end function initialize_symtridiag

   pure module function construct_symtridiag(dv, ev) result(A)
      !! Utility function to create a `SymTridiagonal` matrix from rank-1 arrays.
      real(dp), intent(in) :: dv(:), ev(:)
      !! Diagonal elements of the matrix.
      type(SymTridiagonal) :: A
      !! Corresponding symmetric tridiagonal matrix.
      A%n = size(dv); A%dv = dv; A%ev = ev
      return
   end function construct_symtridiag

   pure module function construct_constant_symtridiag(d, e, n) result(A)
      !! Utility function to create a `SymTridiagonal` matrix with constant diagonal elements.
      real(dp), intent(in) :: d, e
      !! Constant diagonal elements of the matrix.
      integer(ilp), intent(in) :: n
      !! Dimension of the matrix.
      type(SymTridiagonal) :: A
      !! Corresponding symmetric tridiagonal matrix.
      integer i
      A%n = n; A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n-1)]
      return
   end function construct_constant_symtridiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function symtridiag_spmv(A, x) result(y)
      !! Utility function to compute the matrix-vector product where \( A \) is of
      !! `SymTridiagonal` type and `x` and `y` are both rank-1 arrays.
      type(SymTridiagonal), intent(in) :: A
      !! Input matrix.
      real(dp), intent(in) :: x(:)
      !! Input vector.
      real(dp) :: y(size(x))
      !! Output vector.
      integer :: i, n
      n = size(x); y(1) = A%dv(1)*x(1) + A%ev(1)*x(2)
      do concurrent (i=2:n-1)
         y(i) = A%ev(i-1)*x(i-1) + A%dv(i)*x(i) + A%ev(i)*x(i+1)
      enddo
      y(n) = A%dv(n)*x(n) + A%ev(n-1)*x(n-1)
      return
   end function symtridiag_spmv

   pure module function symtridiag_multi_spmv(A, X) result(Y)
      !! Utility function to compute matrix-matrix product where \( A \) is of
      !! `SymTridiagonal` type and `X` and `Y` are both rank-2 arrays.
      type(SymTridiagonal), intent(in) :: A
      !! Input matrix.
      real(dp), intent(in) :: X(:, :)
      !! Input matrix (rank-2 array).
      real(dp) :: Y(size(X, 1), size(X, 2))
      !! Output matrix (rank-2 array).
      integer(ilp) :: i
      do concurrent(i=1:size(X,2))
         Y(:, i) = symtridiag_spmv(A, X(:, i))
      enddo
   end function symtridiag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   pure module function symtridiag_solve(A, b) result(x)
      !! Utility function to solve a linear system \( Ax = b \) where \( A \) is of
      !! `SymTridiagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
      !! with the same type and dimension as `b`.
      type(SymTridiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: b(:)
      !! Right-hand side vector.
      real(dp) :: x(size(b))
      !! Solution vector.
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
      !! Utility function to solve a linear system with multiple right-hand sides where \( A \)
      !! is of `SymTridiagonal` type and `B` a rank-2 array. The solution `X` is also a rank-2
      !! array with the same type and dimension as `B`.
      type(SymTridiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: B(:, :)
      !! Right-hand side vectors.
      real(dp) :: X(size(B, 1), size(B, 2))
      !! Solution vectors.
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
      !! Utility function to convert a `SymTridiagonal` matrix to a regular rank-2 array.
      type(SymTridiagonal), intent(in) :: A
      !! Input symmetric tridiagonal matrix.
      real(dp) :: B(A%n, A%n)
      !! Output dense rank-2 array.
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
      !! Utility function to compute the transpose of a `SymTridiagonal` matrix.
      !! The output matrix is also of `SymTridiagonal` type.
      type(SymTridiagonal), intent(in) :: A
      !! Input symmetric tridiagonal matrix.
      type(SymTridiagonal) :: B
      !! Transpose of the original matrix.
      B = A
      return
   end function symtridiag_transpose

   pure module function symtridiag_shape(A) result(shape)
      type(SymTridiagonal), intent(in) :: A
      !! Input matrix.
      integer(ilp) :: shape(2)
      !! Shape of the matrix.
      shape(2) = A%n
      return
   end function symtridiag_shape

end submodule SymmetricTridiagonal
