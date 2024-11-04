submodule(SpecialMatrices_Tridiagonal) DiagonalMatrices
   use stdlib_linalg, only: eye, diag
   use stdlib_sorting, only: sort, sort_index
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_VALUE_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   module procedure initialize_diag
   ! Utility function to construct a `Diagonal` matrix filled with zeros.
   A%n = n; allocate (A%dv(n)); A%dv = 0.0_dp
   end procedure initialize_diag

   module procedure construct_diag
   ! Utility function to construct a `Diagonal` matrix from a rank-1 array.
   A%n = size(dv); A%dv = dv
   end procedure construct_diag

   module procedure construct_constant_diag
   ! Utility function to construct a `Diagonal` matrix with constant diagonal element.
   integer(ilp) :: i
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
   integer(ilp) :: i
   allocate (y, mold=x)
   do concurrent(i=1:size(x))
      y(i) = A%dv(i)*x(i)
   end do
   end procedure diag_spmv

   module procedure diag_multi_spmv
   ! Utility function to compute the matrix-vector product \( y = A x \) where \( A \)
   ! is of `Diagonal` type and `X` and `Y` are both rank-2 arrays.
   integer(ilp) :: i, j
   allocate (Y, mold=X)
   do concurrent(i=1:size(X, 1), j=1:size(X, 2))
      Y(i, j) = A%dv(i)*X(i, j)
   end do
   end procedure diag_multi_spmv

   module procedure diag_spmv_ip
   integer(ilp) :: i
   do concurrent(i=1:size(x))
      y(i) = A%dv(i)*x(i)
   end do
   return
   end procedure diag_spmv_ip

   module procedure diag_multi_spmv_ip
   integer(ilp) :: i, j
   do concurrent(i=1:size(X, 1), j=1:size(X, 2))
      Y(i, j) = A%dv(i)*X(i, j)
   end do
   return
   end procedure diag_multi_spmv_ip

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   module procedure diag_solve
   ! Utility function to solve the linear system \( A x = b \) where \( A \) is of
   ! `Diagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
   ! with the same type and dimension of `b`.
   integer(ilp) :: i, n
   allocate (x, mold=b)
   do concurrent(i=1:A%n)
      x(i) = b(i)/A%dv(i)
   end do
   end procedure diag_solve

   module procedure diag_multi_solve
   ! Utility function to solve a linear system with multiple right-hande sides where \( A \)
   ! is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a rank-2 array
   ! with the same type and dimensions as `B`.
   integer(ilp) :: i, j
   real(dp) :: inv_dv(A%n)
   allocate (X, mold=B); inv_dv = 1.0_dp/A%dv
   do concurrent(i=1:A%n, j=1:size(B, 2))
      X(i, j) = B(i, j)*inv_dv(i)
   end do
   end procedure diag_multi_solve

   module procedure diag_solve_ip
   ! Utility function to solve *in-place* a linear system with multiple right-hande sides where
   ! \( A \) is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a rank-2 array
   ! with the same type and dimensions as `B`.
   integer(ilp) :: i
   do concurrent(i=1:A%n)
      x(i) = b(i)/A%dv(i)
   end do
   return
   end procedure diag_solve_ip

   module procedure diag_multi_solve_ip
   ! Utility function to solve a linear system with multiple right-hande sides where \( A \)
   ! is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a rank-2 array
   ! with the same type and dimensions as `B`.
   integer(ilp) :: i, j
   real(dp) :: inv_dv(A%n)
   inv_dv = 1.0_dp/A%dv
   do concurrent(i=1:A%n, j=1:size(B, 2))
      X(i, j) = inv_dv(i)*B(i, j)
   end do
   return
   end procedure diag_multi_solve_ip

   module procedure diag_det
   ! Utility function to compute the determinant of a `Diagonal` matrix \( A \).
   determinant = product(A%dv)
   end procedure diag_det

   module procedure diag_trace
   ! Utility function to compute the trace of a `Diagonal` matrix \( A \).
   tr = sum(A%dv)
   end procedure diag_trace

   module procedure diag_inv
   ! Utility function to compute the inverse of a `Diagonal` matrix \( A \).
   B = dense(Diagonal(1.0_dp/A%dv))
   end procedure diag_inv

   module procedure diag_svdvals
   ! Utility function to compute the singular values of a `Diagonal` matrix \( A \).
   s = abs(A%dv); call sort(s, reverse=.true.)
   end procedure diag_svdvals

   module procedure diag_svd
   ! Utility function to compute the singular values of a `Diagonal` matrix \( A \).
   integer(ilp) :: i, index(A%n)
   ! Compute the singular values and sort them in decreasing order.
   s = abs(A%dv); call sort_index(s, index, reverse=.true.)
   ! Compute the left singular vectors.
   u = eye(A%n); u = u(:, index)
   ! Compute the right singular vectors.
   vt = eye(A%n); 
   do concurrent(i=1:A%n, A%dv(i) < 0.0_dp)
      vt(i, i) = -1.0_dp
   end do
   vt = vt(index, :)
   end procedure diag_svd

   module procedure diag_eigvalsh
   ! Utility function to compute the eigenvalues of a real `Diagonal` matrix \( A \).
   lambda = A%dv; call sort(lambda)
   end procedure diag_eigvalsh

   module procedure diag_eigh
   ! Utility function to compute the eigenvalues and eigenvectors of a real `Diagonal`
   ! matrix \( A \).
   integer(ilp) :: index(A%n)
   lambda = A%dv; call sort_index(lambda, index)
   if (present(vectors)) then
      vectors = eye(A%n); vectors = vectors(:, index)
   end if
   end procedure diag_eigh

   !------------------------------------
   !-----     Utility function     -----
   !------------------------------------

   module procedure diag_to_dense
   ! Utility function to convert a `Diagonal` matrix to a regular rank-2 array.
   integer(ilp) :: i
   allocate (B(A%n, A%n)); B = 0.0_dp
   do concurrent(i=1:A%n)
      B(i, i) = A%dv(i)
   end do
   end procedure diag_to_dense

   module procedure diag_transpose
   ! Utility function to compute the transpose of a `Diagonal` matrix.
   B = A
   end procedure diag_transpose

   module procedure diag_shape
   ! Utility function to get the shape of a `Diagonal` matrix.
   shape = A%n
   end procedure diag_shape

   module procedure diag_size
   ! Utility function to get the size of a `Diagonal` matrix along a given dimension.
   arr_size = A%n
   end procedure diag_size

end submodule DiagonalMatrices
