submodule(specialmatrices_circulant) circulant_determinant
   implicit none(type, external)
contains
   module procedure det_rdp
   ! Compute det(A).
   d = real(product(A%c_hat))
   end procedure
end submodule
