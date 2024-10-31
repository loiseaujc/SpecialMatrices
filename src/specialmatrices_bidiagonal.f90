submodule(SpecialMatrices_Tridiagonal) BidiagonalMatrices
   use stdlib_optval, only: optval
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

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   !------------------------------------
   !-----     Utility function     -----
   !------------------------------------

end submodule BidiagonalMatrices
