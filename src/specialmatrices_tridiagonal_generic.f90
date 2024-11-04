submodule(SpecialMatrices_Tridiagonal) GenericTridiagonal
   use stdlib_linalg, only: diag
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   module procedure initialize_tridiag
      !! Utility procedure to construct a `Tridiagonal` matrix filled with zeros.
   A%n = n; allocate (A%d(n), A%du(n - 1), A%dl(n - 1))
   A%d = 0.0_dp; A%du = 0.0_dp; A%dl = 0.0_dp
   return
   end procedure initialize_tridiag

   module procedure construct_tridiag
      !! Utility procedure to construct a `Tridiagonal` matrix from rank-1 arrays.
   A%n = size(d); A%dl = dl; A%d = d; A%du = du
   end procedure construct_tridiag

   module procedure construct_constant_tridiag
      !! Utility procedure to construct a `Tridiagonal` matrix with constant diagonals.
   integer(ilp) :: i
   A%n = n; A%dl = [(l, i=1, n - 1)]; A%d = [(d, i=1, n)]; A%du = [(u, i=1, n - 1)]
   end procedure construct_constant_tridiag

   module procedure construct_dense_to_tridiag
      !! Utility procedure to construct a `Tridiagonal` matrix from a rank-2 array.
   B = Tridiagonal(diag(A, -1), diag(A), diag(A, 1))
   end procedure construct_dense_to_tridiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   module procedure tridiag_spmv
      !! Utility procedure to compute the matrix-vector product \( y = Ax \) where \( A \)
      !! is of `Tridiagonal` type and `x` and `y` are both rank-1 arrays.
   integer(ilp) :: i, n
   n = size(x); allocate (y, mold=x)
   y(1) = A%d(1)*x(1) + A%du(1)*x(2)
   do concurrent(i=2:n - 1)
      y(i) = A%dl(i - 1)*x(i - 1) + A%d(i)*x(i) + A%du(i)*x(i + 1)
   end do
   y(n) = A%d(n)*x(n) + A%dl(n - 1)*x(n - 1)
   end procedure tridiag_spmv

   module procedure tridiag_multi_spmv
      !! Utility procedure to compute the matrix-matrix product \( Y = AX \) where \( A \)
      !!  is of `Tridiagonal` type and `X` and `Y` are both rank-2 arrays.
   integer(ilp) :: i
   allocate (Y, mold=X)
   do concurrent(i=1:size(X, 2))
      Y(:, i) = tridiag_spmv(A, X(:, i))
   end do
   end procedure tridiag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   module procedure tridiag_solve
      !! Utility procedure to solve the linear system \( Ax = b \) where \( A \) is of
      !! `Tridiagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
      !! with the same type and dimension as `b`.
   integer(ilp) :: i, n, nrhs, info
   real(dp) :: dl(A%n - 1), d(A%n), du(A%n - 1), b_(A%n, 1)
   ! Initialize arrays.
   n = A%n; dl = A%dl; d = A%d; du = A%du; b_(:, 1) = b; nrhs = 1
   ! Solve the system.
   call gtsv(n, nrhs, dl, d, du, b_, n, info)
   ! Return the solution.
   x = b_(:, 1)
   end procedure tridiag_solve

   module procedure tridiag_multi_solve
      !! Utility procedure to solve a linear system with multiple right-hand sides where
      !! \( A \) is of `Tridiagonal` type and `B` a rank-2 array. The solution `X` is also
      !! a rank-2 array with the same type and dimensions as `B`.
   integer(ilp) :: i, n, nrhs, info
   real(dp) :: dl(A%n - 1), d(A%n), du(A%n - 1)
   ! Initialize arrays.
   n = A%n; dl = A%dl; d = A%d; du = A%du; nrhs = size(B, 2); X = B
   ! Solve the systems.
   call gtsv(n, nrhs, dl, d, du, X, n, info)
   end procedure tridiag_multi_solve

   !-------------------------------------
   !-----     Utility procedures     -----
   !-------------------------------------

   module procedure tridiag_shape
   shape = A%n
   end procedure tridiag_shape

   module procedure tridiag_to_dense
      !! Utility procedure to convert a `Tridiagonal` matrix to a regular rank-2 array.
   integer(ilp) :: i, n
   n = A%n; B = 0.0_dp
   B(1, 1) = A%d(1); B(1, 2) = A%du(1)
   do concurrent(i=2:n - 1)
      B(i, i - 1) = A%dl(i - 1)
      B(i, i) = A%d(i)
      B(i, i + 1) = A%du(i)
   end do
   B(n, n - 1) = A%dl(n - 1); B(n, n) = A%d(n)
   end procedure tridiag_to_dense

   module procedure tridiag_transpose
      !! Utility procedure to compute the transpose of a `Tridiagonal` matrix.
      !! The output matrix is also of `Tridiagonal` type.
   B = A; B%dl = A%du; B%du = A%dl
   end procedure tridiag_transpose

end submodule GenericTridiagonal
