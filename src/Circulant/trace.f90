submodule(specialmatrices_circulant) circulant_trace
   implicit none(type, external)
contains
   module procedure trace_rdp
   tr = A%c(1) * A%n
   end procedure
end submodule
