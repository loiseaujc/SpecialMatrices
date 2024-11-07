submodule(specialmatrices_strang) strang_determinant
   implicit none(type, external)
contains
   module procedure det_rdp
   d = A%n + 1
   end procedure
end submodule
