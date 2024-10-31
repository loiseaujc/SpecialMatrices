submodule(SpecialMatrices_Tridiagonal) BidiagonalMatrices
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_bidiag(n, which) result(A)
      !! Utility function to construct a `Bidiagonal` matrix filled with zeros.
      integer(ilp), intent(in) :: n
      !! Dimension of the matrix.
      character(len=1), optional, intent(in) :: which
      !! Whether `A` has a sub- or super-diagonal.
      type(Bidiagonal) :: A
      !! Corresponding bidiagonal matrix.
      A%n = n; A%which = optval(which, "L")
      allocate(A%dv(n), A%ev(n-1)); A%dv = 0.0_dp; A%ev = 0.0_dp
      return
   end function initialize_bidiag

   pure module function construct_bidiag(dv, ev, which) result(A)
      !! Utility function to construct a `Bidiagonal` matrix from rank-1 arrays.
      real(dp), intent(in) :: dv(:), ev(:)
      !! Diagonal elements of the matrix.
      character(len=1), optional, intent(in) :: which
      !! Whether `A` has a sub- or super-diagonal.
      type(Bidiagonal) :: A
      !! Corresponding bidiagonal matrix.
      A%n = size(dv); A%dv = dv ; A%ev = ev ; A%which = optval(which, "L")
      return
   end function construct_bidiag

   pure module function construct_constant_bidiag(d, e, n, which) result(A)
      !! Utility function to construct a `Bidiagonal` matrix with constant elements.
      real(dp), intent(in) :: d, e
      !! Constant diagonal elements.
      integer(ilp), intent(in) :: n
      !! Dimension of the matrix.
      character(len=1), optional, intent(in) :: which
      !! Whether `A` has a sub- or super-diagonal.
      type(Bidiagonal) :: A
      !! Corresponding bidiagonal matrix.
      integer :: i
      A%n = n; A%which = optval(which, "L")
      A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n-1)]
      return
   end function construct_constant_bidiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   pure module function bidiag_spmv(A, x) result(y)
      !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
      !! is of type `Bidiagonal` and `x` and `y` are both rank-1 arrays.
      type(Bidiagonal), intent(in) :: A
      !! Input matrix.
      real(dp), intent(in) :: x(:)
      !! Input vector.
      real(dp) :: y(size(x))
      !! Output vector.
      integer(ilp) :: i, n
      n = size(x)
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
      return
   end function bidiag_spmv

   pure module function bidiag_multi_spmv(A, X) result(Y)
      !! Utility function to compute the matrix-matrix product \( Y = AX \) where \( A \)
      !! is of type `Bidiagonal` and `X` and `Y` are both rank-2 arrays.
      type(Bidiagonal), intent(in) :: A
      !! Input matrix.
      real(dp), intent(in) :: X(:, :)
      !! Input matrix (rank-2 array).
      real(dp) :: Y(size(X, 1), size(X, 2))
      !! Output matrix (rank-2 array).
      integer(ilp) :: i, n
      n = size(x)
      do concurrent (i=1:size(X, 2))
         Y(:, i) = bidiag_spmv(A, X(:, i))
      enddo
      return
   end function bidiag_multi_spmv

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   pure module function bidiag_solve(A, b) result(x)
      !! Utility function to solve the linear system \( Ax = b \) where \( A \) is a
      !! `Bidiagonal` matrix and `b` a rank-1 array. The solution `x` is also a rank-1
      !! array with the same type and dimesnion as `b`.
      type(Bidiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: b(:)
      !! Right-hand side vector.
      real(dp) :: x(size(b))
      !! Solution vector.
      real(dp) :: dl(A%n-1), dv(A%n), du(A%n-1), b_(A%n, 1)
      integer :: n, nrhs, info

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
      
      return
   end function bidiag_solve

   pure module function bidiag_multi_solve(A, B) result(X)
      !! Utility function to solve the linear system with multiple right-hand side where 
      !!  \( A \) is a `Bidiagonal` matrix and `B` a rank-2 array. The solution `x` is also
      !! a rank-2 array with the same type and dimesnion as `B`.
      type(Bidiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: B(:, :)
      !! Right-hand side vector.
      real(dp) :: X(size(B, 1), size(B, 2))
      !! Solution vector.
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

      return
   end function bidiag_multi_solve

   !------------------------------------
   !-----     Utility function     -----
   !------------------------------------

   pure module function bidiag_shape(A) result(shape)
      !! Utility function to return the shape of a `Bidiagonal` matrix.
      type(Bidiagonal), intent(in) ::A
      !! Input matrix.
      integer(ilp) :: shape(2)
      !! Shape of the matrix.
      shape = A%n
   end function bidiag_shape

   pure module function bidiag_transpose(A) result(B)
      !! Utility function to compute the transpose of a `Bidiagonal` matrix.
      type(Bidiagonal), intent(in) :: A
      !! Input bidiagonal matrix.
      type(Bidiagonal) :: B
      !! Transpose of the original matrix.
      B = A
      if (A%which == "L") then
         B%which = "U"
      else
         B%which = "L"
      endif
      return
   end function bidiag_transpose

   pure module function bidiag_to_dense(A) result(B)
      !! Utility function to convert a `Bidiagonal` matrix to a regular rank-2 array.
      type(Bidiagonal), intent(in) :: A
      !! Input bidiagonal matrix.
      real(dp) :: B(A%n, A%n)
      !! Output dense rank-2 array.
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
      return
   end function bidiag_to_dense

end submodule BidiagonalMatrices
