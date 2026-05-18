submodule(specialmatrices_hankel) hankel_constructors
   implicit none(type, external)
contains
   module procedure construct
   !> Dimension of the matrix.
   A%m = size(vc) ; A%n = size(vr)
   !> First column/Last row vectors.
   A%vc = vc ; A%vr = vr
   !> Ensure vc[end] and vr[1] are the same.
   A%vr(A%m) = A%vc(1)
   end procedure construct
end submodule hankel_constructors
