submodule(SpecialMatrices_Tridiagonal) DiagonalMatrices
   use stdlib_linalg, only: eye, diag
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_diag(n) result(A)
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Output matrix.
      type(Diagonal) :: A
      A%n = n; allocate(A%dv(n)); A%dv = 0.0_wp
      return
   end function initialize_diag

   pure module function construct_diag(dv) result(A)
      ! Diagonal elements.
      real(wp), intent(in) :: dv(:)
      ! Output matrix.
      type(Diagonal) :: A
      A%n = size(dv); A%dv = dv
      return
   end function construct_diag

   pure module function construct_constant_diag(d, n) result(A)
      ! Diagonal element.
      real(wp), intent(in) :: d
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Output matrix.
      type(Diagonal) :: A
      integer :: i
      A%n = n; A%dv = [(d, i=1, n)]
      return
   end function construct_constant_diag

   module function construct_dense_to_diag(A) result(B)
      ! Input dense matrix.
      real(wp), intent(in) :: A(:, :)
      ! Output diagonal matrix.
      type(Diagonal) :: B
      B = Diagonal(diag(A))
      return
   end function construct_dense_to_diag

   !------------------------------------------------------------
   !-----     Matrix-vector and matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function diag_spmv(A, x) result(y)
      ! Input matrix.
      type(Diagonal), intent(in) :: A
      ! Input vector.
      real(wp), intent(in) :: x(:)
      ! Output vector.
      real(wp) :: y(size(x))
      integer :: i
      do concurrent(i=1:size(x))
         y(i) = A%dv(i) * x(i)
      enddo
      return
   end function diag_spmv

   pure module function diag_multi_spmv(A, X) result(Y)
      ! Input matrix.
      type(Diagonal), intent(in) :: A
      ! Iput vectors.
      real(wp), intent(in) :: X(:, :)
      ! Output vectors.
      real(wp) :: Y(size(X, 1), size(X, 2))
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
      ! Coefficient matrix.
      type(Diagonal), intent(in) :: A
      ! Right-hande side vector.
      real(wp), intent(in) :: b(:)
      ! Solution vector.
      real(wp) :: x(size(b))
      integer :: i
      do concurrent(i=1:A%n)
         x(i) = b(i) / A%dv(i)
      enddo
      return
   end function diag_solve

   pure module function diag_multi_solve(A, B) result(X)
      ! Coefficient matrix.
      type(Diagonal), intent(in) :: A
      ! Right-hand side vectors.
      real(wp), intent(in) :: B(:, :)
      ! Solution vectors.
      real(wp) :: X(size(B, 1), size(B, 2))
      integer(int32) :: i, j
      real(wp) :: dv(A%n)
      dv = 1.0_wp / A%dv
      do concurrent(i=1:A%n, j=1:size(B, 2))
         X(i, j) = B(i, j) / dv(i)
      enddo
      return
   end function diag_multi_solve

   ! module subroutine diag_eig(A, lambda, vectors)
   !    ! Input matrix.
   !    type(Diagonal), intent(in) :: A
   !    ! Eigenvalues.
   !    real(wp), intent(out) :: lambda(A%n)
   !    ! Eigenvectors.
   !    real(wp), intent(out) :: vectors(A%n, A%n)
   !
   !    lambda = A%dv ; vectors = eye(A%n)
   !    return
   ! end subroutine diag_eig

   !------------------------------------
   !-----     Utility function     -----
   !------------------------------------

   pure module function diag_to_dense(A) result(B)
      ! Input Diagonal matrix.
      type(Diagonal), intent(in) :: A
      ! Output dense matrix.
      real(wp) :: B(A%n, A%n)
      integer :: i
      B = 0.0_wp
      do concurrent(i=1:A%n)
         B(i, i) = A%dv(i)
      enddo
   end function diag_to_dense

   pure module function diag_transpose(A) result(B)
      ! Input Diagonal matrix.
      type(Diagonal), intent(in) :: A
      ! Output matrix.
      type(Diagonal) :: B
      B = A
   end function diag_transpose
   
end submodule DiagonalMatrices
