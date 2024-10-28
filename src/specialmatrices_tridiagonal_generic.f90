submodule(SpecialMatrices_Tridiagonal) GenericTridiagonal
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_tridiag(n) result(A)
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Output matrix.
      type(Tridiagonal) :: A
      A%n = n
      allocate (A%d(n), A%du(n - 1), A%dl(n - 1))
      A%d = 0.0_wp; A%du = 0.0_wp; A%dl = 0.0_wp
      return
   end function initialize_tridiag

   pure module function construct_tridiag(dl, d, du) result(A)
      ! Diagonals of the matrix.
      real(wp), intent(in) :: dl(:), d(:), du(:)
      ! Output matrix.
      type(Tridiagonal) :: A
      A%n = size(d); A%dl = dl; A%d = d; A%du = du
      return
   end function construct_tridiag

   pure module function construct_constant_tridiag(l, d, u, n) result(A)
      ! Constant diagonals.
      real(wp), intent(in) :: l, d, u
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Output matrix.
      type(Tridiagonal) :: A
      integer(int32) :: i
      A%n = n; A%dl = [(l, i=1, n - 1)]; A%d = [(d, i=1, n)]; A%du = [(u, i=1, n - 1)]
   end function construct_constant_tridiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function tridiag_spmv(A, x) result(y)
      ! Input matrix.
      type(Tridiagonal), intent(in) :: A
      ! Input vector.
      real(wp), intent(in) :: x(:)
      ! Output vector.
      real(wp) :: y(size(x))
      integer :: i, n
      n = size(x); y(1) = A%d(1)*x(1) + A%du(1)*x(2)
      do concurrent(i=2:n - 1)
         y(i) = A%dl(i)*x(i - 1) + A%d(i)*x(i) + A%du(i)*x(i + 1)
      end do
      y(n) = A%d(n)*x(n) + A%dl(n - 1)*x(n - 1)
      return
   end function tridiag_spmv

   pure module function tridiag_multi_spmv(A, X) result(Y)
      ! Input matrix.
      type(Tridiagonal), intent(in) :: A
      ! Inputs vectors.
      real(wp), intent(in) :: X(:, :)
      ! Output vectors.
      real(wp) :: Y(size(X, 1), size(X, 2))
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
      ! Coefficient matrix.
      type(Tridiagonal), intent(in) :: A
      ! Right-hand side vector.
      real(wp), intent(in) :: b(:)
      ! Solution vector.
      real(wp) :: x(size(b))
      integer(int32) :: i, n
      real(wp) :: w, dl(A%n - 1), d(A%n), du(A%n - 1), b_(A%n)

      ! Initialize arrays.
      n = A%n; dl = A%dl; d = A%d; du = A%du; b_ = b

      ! Update matrix coefficients.
      do i = 2, n
         w = dl(i - 1)/d(i - 1)
         d(i) = d(i) - w*du(i - 1)
         b_(i) = b_(i) - w*b_(i - 1)
      end do

      ! Backward substitution.
      x(n) = b_(n)/d(n)
      do i = n - 1, 1
         x(i) = (b_(i) - du(i)*x(i + 1))/d(i)
      end do

      return
   end function tridiag_solve

   pure module function tridiag_multi_solve(A, B) result(X)
      ! Coefficient matrix.
      type(Tridiagonal), intent(in) :: A
      ! Right-hand side vectors.
      real(wp), intent(in) :: B(:, :)
      ! Solution vectors.
      real(wp) :: X(size(B, 1), size(B, 2))
      integer :: i
      do concurrent(i=1:size(B, 2))
         X(:, i) = tridiag_solve(A, B(:, i))
      end do
   end function tridiag_multi_solve

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   pure module function tridiag_to_dense(A) result(B)
      ! Input tridiagonal matrix.
      type(Tridiagonal), intent(in) :: A
      ! Output dense matrix.
      real(wp) :: B(A%n, A%n)
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
      ! Input tridiagonal matrix.
      type(Tridiagonal), intent(in) :: A
      ! Transposed matrix.
      type(Tridiagonal) :: B
      B = A; B%dl = A%du; B%du = A%dl
   end function tridiag_transpose

end submodule GenericTridiagonal
