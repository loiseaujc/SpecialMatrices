submodule(specialmatrices_circulant) circulant_inverse
   implicit none(type, external)
contains
   module procedure inv_rdp
   B = circulant(real(ifft(1.0_dp/A%c_hat)) / A%n)
   end procedure inv_rdp
end submodule circulant_inverse
