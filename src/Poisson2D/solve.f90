submodule(specialmatrices_poisson2D) poisson2D_solve
   implicit none(type, external)
contains
   module procedure solve_single_rhs
   end procedure

   module procedure solve_multi_rhs
   end procedure
end submodule
