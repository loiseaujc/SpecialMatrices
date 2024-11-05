submodule(specialmatrices_symtridiagonal) symtridiagonal_trace
   implicit none(type, external)
contains
   module procedure trace_rdp
   tr = sum(A%dv)
   end procedure
end submodule
