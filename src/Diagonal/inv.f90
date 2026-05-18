submodule(specialmatrices_diagonal) diagonal_inverse
   implicit none(type, external)
contains
   module procedure inv_rdp
   B = Diagonal(1.0_dp/A%dv)
   end procedure inv_rdp
end submodule diagonal_inverse
