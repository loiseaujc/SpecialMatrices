submodule(SpecialMatrices_Tridiagonal) DiagonalMatrices
   use stdlib_linalg, only: eye, diag
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_diag(n) result(A)
      !! Utility function to construct a `Diagonal` matrix filled with zeros.
      integer(int32), intent(in) :: n
      !! Dimension of the matrix.
      type(Diagonal) :: A
      !! Corresponding diagonal matrix.
      A%n = n; allocate(A%dv(n)); A%dv = 0.0_wp
      return
   end function initialize_diag

   pure module function construct_diag(dv) result(A)
      !! Utility function to construct a `Diagonal` matrix from a rank-1 array.
      real(wp), intent(in) :: dv(:)
      !! Diagonal element of the matrix.
      type(Diagonal) :: A
      !! Corresponding diagonal matrix.
      A%n = size(dv); A%dv = dv
      return
   end function construct_diag

   pure module function construct_constant_diag(d, n) result(A)
      !! Utility function to construct a `Diagonal` matrix with constant diagonal element.
      real(wp), intent(in) :: d
      !! Constant diagonal element of the matrix.
      integer(int32), intent(in) :: n
      !! Dimension of the matrix.
      type(Diagonal) :: A
      !! Corresponding diagonal matrix.
      integer :: i
      A%n = n; A%dv = [(d, i=1, n)]
      return
   end function construct_constant_diag

   module function construct_dense_to_diag(A) result(B)
      !! Utility function to construct a `Diagonal` matrix from a rank-2 array.
      real(wp), intent(in) :: A(:, :)
      !! Dense [n x n] matrix from which to construct the `Diagonal` one.
      type(Diagonal) :: B
      !! Corresponding diagonal matrix.
      B = Diagonal(diag(A))
      return
   end function construct_dense_to_diag

   !------------------------------------------------------------
   !-----     Matrix-vector and matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function diag_spmv(A, x) result(y)
      !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
      !! is of `Diagonal` type and `x` and `y` are both rank-1 arrays.
      type(Diagonal), intent(in) :: A
      !! Input matrix.
      real(wp), intent(in) :: x(:)
      !! Input vector.
      real(wp) :: y(size(x))
      !! Output vector.
      integer :: i
      do concurrent(i=1:size(x))
         y(i) = A%dv(i) * x(i)
      enddo
      return
   end function diag_spmv

   pure module function diag_multi_spmv(A, X) result(Y)
      !! Utility function to compute the matrix-vector product \( y = A x \) where \( A \)
      !! is of `Diagonal` type and `X` and `Y` are both rank-2 arrays.
      type(Diagonal), intent(in) :: A
      !! Input matrix.
      real(wp), intent(in) :: X(:, :)
      !! Input matrix (rank-2 array).
      real(wp) :: Y(size(X, 1), size(X, 2))
      !! Output matrix (rank-2 array).
      integer :: i, j
      do concurrent(i=1:size(X, 1), j=1:size(X, 2))
         Y(i, j) = A%dv(i) * X(i, j)
      enddo
      return
   end function diag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   pure module function diag_solve(A, b) result(x)
      !! Utility function to solve the linear system \( A x = b \) where \( A \) is of
      !! `Diagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
      !! with the same type and dimension of `b`.
      type(Diagonal), intent(in) :: A
      !! Coefficient matrix.
      real(wp), intent(in) :: b(:)
      !! Right-hande side vector.
      real(wp) :: x(size(b))
      !! Solution vector.
      integer :: i
      do concurrent(i=1:A%n)
         x(i) = b(i) / A%dv(i)
      enddo
      return
   end function diag_solve

   pure module function diag_multi_solve(A, B) result(X)
      !! Utility function to solve a linear system with multiple right-hande sides where \( A \)
      !! is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a rank-2 array
      !! with the same type and dimensions as `B`.
      type(Diagonal), intent(in) :: A
      !! Coefficient matrix.
      real(wp), intent(in) :: B(:, :)
      !! Right-hand side vectors.
      real(wp) :: X(size(B, 1), size(B, 2))
      !! Solution vectors.
      integer(int32) :: i, j
      real(wp) :: dv(A%n)
      dv = 1.0_wp / A%dv
      do concurrent(i=1:A%n, j=1:size(B, 2))
         X(i, j) = B(i, j) / dv(i)
      enddo
      return
   end function diag_multi_solve

   pure module subroutine diag_eig(A, lambda, vectors)
      type(Diagonal), intent(in) :: A
      real(wp), intent(out) :: lambda(A%n)
      real(wp), intent(out) :: vectors(A%n, A%n)

      lambda = A%dv ; vectors = eye(A%n)
      return
   end subroutine diag_eig

   !------------------------------------
   !-----     Utility function     -----
   !------------------------------------

   pure module function diag_to_dense(A) result(B)
      !! Utility function to convert a `Diagonal` matrix to a regular rank-2 array.
      type(Diagonal), intent(in) :: A
      !! Input diagonal matrix.
      real(wp) :: B(A%n, A%n)
      !! Output dense rank-2 array.
      integer :: i
      B = 0.0_wp
      do concurrent(i=1:A%n)
         B(i, i) = A%dv(i)
      enddo
   end function diag_to_dense

   pure module function diag_transpose(A) result(B)
      !! Utility function to compute the transpose of a `Diagonal` matrix.
      type(Diagonal), intent(in) :: A
      !! Input matrix.
      type(Diagonal) :: B
      !! Transposed matrix.
      B = A
   end function diag_transpose
   
end submodule DiagonalMatrices
