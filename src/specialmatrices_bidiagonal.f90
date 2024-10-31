submodule(SpecialMatrices_Tridiagonal) BidiagonalMatrices
   use stdlib_optval, only: optval
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_bidiag(n, which) result(A)
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Upper- or lower-bidiagonal.
      character(len=1), optional, intent(in) :: which
      ! Output matrix.
      type(Bidiagonal) :: A
      A%n = n; A%which = optval(which, "L")
      allocate(A%dv(n), A%ev(n-1)); A%dv = 0.0_wp; A%ev = 0.0_wp
      return
   end function initialize_bidiag

   pure module function construct_bidiag(dv, ev, which) result(A)
      ! Diagonal elements.
      real(wp), intent(in) :: dv(:), ev(:)
      ! Upper- or lower-bidiagonal.
      character(len=1), optional, intent(in) :: which
      ! Output matrix.
      type(Bidiagonal) :: A
      A%n = size(dv); A%dv = dv ; A%ev = ev ; A%which = optval(which, "L")
      return
   end function construct_bidiag

   pure module function construct_constant_bidiag(d, e, n, which) result(A)
      ! Constant diagonal elements.
      real(wp), intent(in) :: d, e
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Upper- or lower-bidiagonal.
      character(len=1), optional, intent(in) :: which
      ! Output matrix.
      type(Bidiagonal) :: A
      integer :: i
      A%n = n; A%which = optval(which, "L")
      A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n-1)]
      return
   end function construct_constant_bidiag
end submodule BidiagonalMatrices
