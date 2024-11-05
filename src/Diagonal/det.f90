submodule(specialmatrices_diagonal) determinant
   implicit none(type, external)
contains
   module procedure det_rdp
   ! Compute det(A).
   d = product(A%dv)
   end procedure
end submodule
