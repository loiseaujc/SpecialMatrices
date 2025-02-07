submodule(specialmatrices_poisson2D) poisson2D_eigh
   use stdlib_constants, only: pi => pi_dp
   use stdlib_sorting, only: sort_index
   implicit none(type, external)
contains
   module procedure eigvalsh_rdp
   end procedure

   module procedure eigh_rdp
   end procedure
end submodule
