submodule(SpecialMatrices_Tridiagonal) DiagonalMatrices
   use stdlib_linalg, only: eye, diag
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_VALUE_ERROR, &
      LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   module procedure initialize_diag
      ! Utility function to construct a `Diagonal` matrix filled with zeros.
      A%n = n; allocate(A%dv(n)); A%dv = 0.0_dp
   end procedure initialize_diag

   module procedure construct_diag
      ! Utility function to construct a `Diagonal` matrix from a rank-1 array.
      A%n = size(dv); A%dv = dv
   end procedure construct_diag

   module procedure construct_constant_diag
      ! Utility function to construct a `Diagonal` matrix with constant diagonal element.
      integer :: i
      A%n = n; A%dv = [(d, i=1, n)]
   end procedure construct_constant_diag

   module procedure construct_dense_to_diag
      ! Utility function to construct a `Diagonal` matrix from a rank-2 array.
      B = Diagonal(diag(A))
   end procedure construct_dense_to_diag

   !------------------------------------------------------------
   !-----     Matrix-vector and matrix-matrix products     -----
   !------------------------------------------------------------

   module procedure diag_spmv
      ! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
      ! is of `Diagonal` type and `x` and `y` are both rank-1 arrays.
      integer :: i
      allocate(y, mold=x)
      do concurrent(i=1:size(x))
         y(i) = A%dv(i) * x(i)
      enddo
   end procedure diag_spmv

   module procedure diag_multi_spmv
      ! Utility function to compute the matrix-vector product \( y = A x \) where \( A \)
      ! is of `Diagonal` type and `X` and `Y` are both rank-2 arrays.
      integer :: i, j
      allocate(Y, mold=X)
      do concurrent(i=1:size(X, 1), j=1:size(X, 2))
         Y(i, j) = A%dv(i) * X(i, j)
      enddo
   end procedure diag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   module procedure diag_solve
      ! Utility function to solve the linear system \( A x = b \) where \( A \) is of
      ! `Diagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
      ! with the same type and dimension of `b`.
      integer :: i, n
      if (.not. allocated(x)) allocate(x, mold=b)
      do concurrent(i=1:A%n)
         x(i) = b(i) / A%dv(i)
      enddo
   end procedure diag_solve

   module procedure diag_multi_solve
      ! Utility function to solve a linear system with multiple right-hande sides where \( A \)
      ! is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a rank-2 array
      ! with the same type and dimensions as `B`.
      integer(ilp) :: i, j
      real(dp) :: inv_dv(A%n)
      if (.not. allocated(X)) allocate(X, mold=B); inv_dv = 1.0_dp / A%dv
      do concurrent(i=1:A%n, j=1:size(B, 2))
         X(i, j) = B(i, j) * inv_dv(i)
      enddo
   end procedure diag_multi_solve

   module procedure diag_eig
      lambda = A%dv ; vectors = eye(A%n)
   end procedure diag_eig

   !------------------------------------
   !-----     Utility function     -----
   !------------------------------------

   module procedure diag_to_dense
      ! Utility function to convert a `Diagonal` matrix to a regular rank-2 array.
      integer :: i
      B = 0.0_dp
      do concurrent(i=1:A%n)
         B(i, i) = A%dv(i)
      enddo
   end procedure diag_to_dense

   module procedure diag_transpose
      ! Utility function to compute the transpose of a `Diagonal` matrix.
      B = A
   end procedure diag_transpose

   module procedure diag_shape
      ! Utility function to get the shape of a `Diagonal` matrix.
      shape = A%n
   end procedure diag_shape
   
end submodule DiagonalMatrices
