submodule(specialmatrices_poisson2D) poisson2D_trace
   implicit none(type, external)
contains
   module procedure trace_rdp
   integer(ilp) :: n
   n = size(A, 1)
   tr = -n * (2/A%dx**2 + 2/A%dy**2)
   end procedure
end submodule
