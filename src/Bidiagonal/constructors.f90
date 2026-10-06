submodule(specialmatrices_bidiagonal) bidiagonal_constructors
   use stdlib_optval, only: optval
   use stdlib_linalg, only: diag
   implicit none(type, external)
contains
   module procedure initialize
   A%n = n; allocate (A%dv(n)); allocate (A%ev(n - 1))
   A%dv = 0.0_dp; A%ev = 0.0_dp; A%which = "L"
   end procedure initialize

   module procedure construct
   A%n = size(dv)
   allocate (A%dv(A%n), source=0.0_dp); A%dv = dv
   allocate (A%ev(A%n - 1), source=0.0_dp); A%ev = ev
   A%which = optval(which, "L")
   end procedure construct

   module procedure construct_constant
   integer(ilp) :: i
   A%n = n; A%which = optval(which, "L")
   A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n - 1)]
   end procedure construct_constant
end submodule bidiagonal_constructors
