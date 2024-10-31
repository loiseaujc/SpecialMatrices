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
end submodule BidiagonalMatrices
