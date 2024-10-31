submodule(SpecialMatrices_Tridiagonal) GenericTridiagonal
   use stdlib_linalg, only: diag
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_tridiag(n) result(A)
      !! Utility function to construct a `Tridiagonal` matrix filled with zeros.
      integer(int32), intent(in) :: n
      !! Dimension of the matrix.
      type(Tridiagonal) :: A
      !! Corresponding tridiagonal matrix.
      A%n = n; allocate (A%d(n), A%du(n - 1), A%dl(n - 1))
      A%d = 0.0_wp; A%du = 0.0_wp; A%dl = 0.0_wp
      return
   end function initialize_tridiag

   pure module function construct_tridiag(dl, d, du) result(A)
      !! Utility function to construct a `Tridiagonal` matrix from rank-1 arrays.
      real(wp), intent(in) :: dl(:), d(:), du(:)
      !! Diagonal elements of the matrix.
      type(Tridiagonal) :: A
      !! Corresponding tridiagonal matrix.
      A%n = size(d); A%dl = dl; A%d = d; A%du = du
      return
   end function construct_tridiag

   pure module function construct_constant_tridiag(l, d, u, n) result(A)
      !! Utility function to construct a `Tridiagonal` matrix with constant diagonals.
      real(wp), intent(in) :: l, d, u
      !! Constant diagonal elements of the matrix.
      integer(int32), intent(in) :: n
      !! Dimension of the matrix.
      type(Tridiagonal) :: A
      !! Corresponding tridiagonal matrix.
      integer(int32) :: i
      A%n = n; A%dl = [(l, i=1, n - 1)]; A%d = [(d, i=1, n)]; A%du = [(u, i=1, n - 1)]
   end function construct_constant_tridiag

   module function construct_dense_to_tridiag(A) result(B)
      !! Utility function to construct a `Tridiagonal` matrix from a rank-2 array.
      real(wp), intent(in) :: A(:, :)
      !! Dense [n x n] matrix from which to construct the `Tridiagonal` one.
      type(Tridiagonal) :: B
      !! Corresponding tridiagonal matrix.
      B = Tridiagonal(diag(A, -1), diag(A), diag(A, 1))
      return
   end function construct_dense_to_tridiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function tridiag_spmv(A, x) result(y)
      !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
      !! is of `Tridiagonal` type and `x` and `y` are both rank-1 arrays.
      type(Tridiagonal), intent(in) :: A
      !! Input matrix.
      real(wp), intent(in) :: x(:)
      !! Input vector.
      real(wp) :: y(size(x))
      !! Output vector.
      integer :: i, n
      n = size(x); y(1) = A%d(1)*x(1) + A%du(1)*x(2)
      do concurrent(i=2:n - 1)
         y(i) = A%dl(i)*x(i - 1) + A%d(i)*x(i) + A%du(i)*x(i + 1)
      end do
      y(n) = A%d(n)*x(n) + A%dl(n - 1)*x(n - 1)
      return
   end function tridiag_spmv

   pure module function tridiag_multi_spmv(A, X) result(Y)
      !! Utility function to compute the matrix-matrix product \( Y = AX \) where \( A \)
      !!  is of `Tridiagonal` type and `X` and `Y` are both rank-2 arrays.
      type(Tridiagonal), intent(in) :: A
      !! Input matrix.
      real(wp), intent(in) :: X(:, :)
      !! Input matrix (rank-2 array).
      real(wp) :: Y(size(X, 1), size(X, 2))
      !! Output matrix (rank-2 array).
      integer(int32) :: i
      do concurrent(i=1:size(X, 2))
         Y(:, i) = tridiag_spmv(A, X(:, i))
      end do
      return
   end function tridiag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   pure module function tridiag_solve(A, b) result(x)
      !! Utility function to solve the linear system \( Ax = b \) where \( A \) is of
      !! `Tridiagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
      !! with the same type and dimension as `b`.
      type(Tridiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(wp), intent(in) :: b(:)
      !! Right-hand side vector.
      real(wp) :: x(size(b))
      !! Solution vector.
      integer :: i, n, nrhs, info
      real(wp) :: dl(A%n - 1), d(A%n), du(A%n - 1), b_(A%n, 1)

      ! Initialize arrays.
      n = A%n; dl = A%dl; d = A%d; du = A%du; b_(:, 1) = b ; nrhs = 1

      ! Solve the system.
      call gtsv(n, nrhs, dl, d, du, b_, n, info)

      ! Return the solution.
      x = b_(:, 1)

      return
   end function tridiag_solve

   pure module function tridiag_multi_solve(A, B) result(X)
      !! Utility function to solve a linear system with multiple right-hand sides where
      !! \( A \) is of `Tridiagonal` type and `B` a rank-2 array. The solution `X` is also
      !! a rank-2 array with the same type and dimensions as `B`.
      type(Tridiagonal), intent(in) :: A
      real(wp), intent(in) :: B(:, :)
      real(wp) :: X(size(B, 1), size(B, 2))
      integer :: i, n, nrhs, info
      real(wp) :: dl(A%n - 1), d(A%n), du(A%n - 1)

      ! Initialize arrays.
      n = A%n; dl = A%dl; d = A%d; du = A%du; nrhs = size(B, 2); X = B

      ! Solve the systems.
      call gtsv(n, nrhs, dl, d, du, X, n, info)

   end function tridiag_multi_solve

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   pure module function tridiag_to_dense(A) result(B)
      !! Utility function to convert a `Tridiagonal` matrix to a regular rank-2 array.
      type(Tridiagonal), intent(in) :: A
      !! Input tridiagonal matrix.
      real(wp) :: B(A%n, A%n)
      !! Output dense rank-2 array.
      integer :: i, n
      n = A%n; B = 0.0_wp
      B(1, 1) = A%d(1); B(1, 2) = A%du(1)
      do concurrent(i=2:n - 1)
         B(i, i - 1) = A%dl(i - 1)
         B(i, i) = A%d(i)
         B(i, i + 1) = A%du(i)
      end do
      B(n, n - 1) = A%dl(n - 1); B(n, n) = A%d(n)
      return
   end function tridiag_to_dense

   pure module function tridiag_transpose(A) result(B)
      !! Utility function to compute the transpose of a `Tridiagonal` matrix.
      !! The output matrix is also of `Tridiagonal` type.
      type(Tridiagonal), intent(in) :: A
      !! Input tridiagonal matrix.
      type(Tridiagonal) :: B
      !! Transpose of the original tridiagonal matrix.
      B = A; B%dl = A%du; B%du = A%dl
   end function tridiag_transpose

end submodule GenericTridiagonal
