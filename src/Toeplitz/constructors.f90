submodule(specialmatrices_toeplitz) toeplitz_constructors
   implicit none(type, external)
contains
   module procedure construct
   !> Dimension of the matrix.
   A%m = size(vc) ; A%n = size(vr)
   !> First column/row vectors.
   A%vc = vc ; A%vr = vr
   !> Ensure vc[1] and vr[1] are the same.
   A%vr(1) = A%vc(1)
   end procedure
end submodule
