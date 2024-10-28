submodule(SpecialMatrices_Tridiagonal) SymmetricTridiagonal
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   pure module function initialize_symtridiag(n) result(A)
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Output matrix.
      type(SymTridiagonal) :: A
      A%n = n; allocate(A%dv(n), A%ev(n-1))
      A%dv = 0.0_wp; A%ev = 0.0_wp
      return
   end function initialize_symtridiag

   pure module function construct_symtridiag(dv, ev) result(A)
      ! Diagonals elements.
      real(wp), intent(in) :: dv(:), ev(:)
      ! Output matrix.
      type(SymTridiagonal) :: A
      A%n = size(dv); A%dv = dv; A%ev = ev
      return
   end function construct_symtridiag

   pure module function construct_constant_symtridiag(d, e, n) result(A)
      ! Diagonal elements.
      real(wp), intent(in) :: d, e
      ! Dimension of the matrix.
      integer(int32), intent(in) :: n
      ! Output matrix.
      type(SymTridiagonal) :: A
      integer i
      A%n = n; A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n-1)]
      return
   end function construct_constant_symtridiag

end submodule SymmetricTridiagonal
