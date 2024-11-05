submodule(specialmatrices_tridiagonal) tridiagonal_constructors
   implicit none(type, external)
contains
   module procedure initialize
   A%n = n; allocate (A%dl(n - 1)); allocate (A%dv(n)); allocate (A%du(n-1))
   A%dl = 0.0_dp; A%dv = 0.0_dp; A%du = 0.0_dp 
   end procedure

   module procedure construct
   integer(ilp) :: n
   n = size(dv)
   A%n = n; A%dl = dl; A%dv = dv; A%du = du
   end procedure

   module procedure construct_constant
   integer(ilp) :: i
   A%n = n
   A%dl = [(dl, i=1, n - 1)]; A%dv = [(dv, i=1, n)]; A%du = [(du, i=1, n-1)]
   end procedure
end submodule
