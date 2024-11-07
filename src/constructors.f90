submodule(specialmatrices_strang) strang_constructors
   implicit none(type, external)
contains
   module procedure initialize
   A%n = n
   end procedure
end submodule
