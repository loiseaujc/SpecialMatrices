submodule(SpecialMatrices_Tridiagonal) BidiagonalMatrices
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   module procedure initialize_bidiag
      ! Utility procedure to construct a `Bidiagonal` matrix filled with zeros.
      A%n = n; A%which = optval(which, "L")
      allocate(A%dv(n), A%ev(n-1)); A%dv = 0.0_dp; A%ev = 0.0_dp
   end procedure initialize_bidiag

   module procedure construct_bidiag
      !! Utility procedure to construct a `Bidiagonal` matrix from rank-1 arrays.
      A%n = size(dv); A%dv = dv ; A%ev = ev ; A%which = optval(which, "L")
   end procedure construct_bidiag

   module procedure construct_constant_bidiag
      !! Utility procedure to construct a `Bidiagonal` matrix with constant elements.
      integer(ilp) :: i
      A%n = n; A%which = optval(which, "L")
      A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n-1)]
   end procedure construct_constant_bidiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   module procedure bidiag_spmv
      integer(ilp) :: i, n
      n = size(x); allocate(y, mold=x)
      if (A%which == "L") then
         y(1) = A%dv(1) * x(1)
         do concurrent (i=2:n)
            y(i) = A%ev(i-1)*x(i-1) + A%dv(i)*x(i)
         enddo
      elseif (A%which == "U") then
         do concurrent (i=1:n-1)
            y(i) = A%dv(i)*x(i) + A%ev(i)*x(i+1)
         enddo
         y(n) = A%dv(n) * x(n)
      endif
   end procedure bidiag_spmv

   module procedure bidiag_multi_spmv
      integer(ilp) :: i, n
      n = size(x); allocate(Y, mold=X)
      do concurrent (i=1:size(X, 2))
         Y(:, i) = bidiag_spmv(A, X(:, i))
      enddo
   end procedure bidiag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   module procedure bidiag_solve
      real(dp) :: dl(A%n-1), dv(A%n), du(A%n-1), b_(A%n, 1)
      integer(ilp) :: n, nrhs, info
      ! Initialize array.
      n = A%n; b_(:, 1) = b; nrhs = 1; dv = A%dv
      ! Dispatch based on upper- or lower-bidiagonal.
      if (A%which == "L") then
         dl = A%ev; du = 0.0_dp
      elseif (A%which == "U") then
         dl = 0.0_dp; du = A%ev
      endif
      ! Solve the system using a tridiagonal solver.
      call gtsv(n, nrhs, dl, dv, du, b_, n, info)
      ! Return the solution.
      x = b_(:, 1)
   end procedure bidiag_solve

   module procedure bidiag_multi_solve
      real(dp) :: dl(A%n-1), dv(A%n), du(A%n-1), b_(A%n, 1)
      integer(ilp) :: n, nrhs, info
      ! Initialize array.
      n = A%n; X = B; nrhs = size(X, 2); dv = A%dv
      ! Dispatch based on upper- or lower-bidiagonal.
      if (A%which == "L") then
         dl = A%ev; du = 0.0_dp
      elseif (A%which == "U") then
         dl = 0.0_dp; du = A%ev
      endif
      ! Solve the system using a tridiagonal solver.
      call gtsv(n, nrhs, dl, dv, du, X, n, info)
   end procedure bidiag_multi_solve

   !------------------------------------
   !-----     Utility procedure     -----
   !------------------------------------

   module procedure bidiag_shape
      shape = A%n
   end procedure bidiag_shape

   module procedure bidiag_transpose
      B = A
      if (A%which == "L") then
         B%which = "U"
      else
         B%which = "L"
      endif
   end procedure bidiag_transpose

   module procedure bidiag_to_dense
      integer(ilp) :: i, n
      n = A%n; B = 0.0_dp
      if (A%which == "L") then
         B(1, 1) = A%dv(1)
         do concurrent(i=2:n)
            B(i, i) = A%dv(i)
            B(i, i-1) = A%ev(i-1)
         enddo
      elseif (A%which == "U") then
         do concurrent(i=1:n-1)
            B(i, i) = A%dv(i)
            B(i, i+1) = A%ev(i)
         enddo
         B(n, n) = A%dv(n)
      endif
   end procedure bidiag_to_dense

end submodule BidiagonalMatrices
