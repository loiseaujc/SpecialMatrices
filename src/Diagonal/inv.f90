submodule(specialmatrices_diagonal) inverse
   implicit none(type, external)
contains
   module procedure inv_rdp
   B = Diagonal(1.0_dp/A%dv)
   end procedure
end submodule
