submodule(specialmatrices_bidiagonal) bidiagonal_determinant
   implicit none(type, external)
contains
   module procedure det_rdp
   d = product(A%dv)
   end procedure det_rdp
end submodule bidiagonal_determinant
