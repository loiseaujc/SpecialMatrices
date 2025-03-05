submodule(specialmatrices_diagonal) diagonal_trace
   implicit none(type, external)
contains
   module procedure trace_rdp
   tr = sum(A%dv)
   end procedure
end submodule
