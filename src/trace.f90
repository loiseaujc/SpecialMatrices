submodule(specialmatrices_strang) strange_trace
   implicit none(type, external)
contains
   module procedure trace_rdp
   tr = A%n*2
   end procedure
end submodule
