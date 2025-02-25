submodule(specialmatrices_poisson2D) poisson2D_constructors
   use stdlib_optval, only: optval
   implicit none(type, external)
contains
   module procedure initialize
   ! Grid spacing.
   real(dp) :: dx, dy
   A%nx = nx; A%ny = ny;
   dx = optval(Lx, 1.0_dp)/(nx+1); dy = optval(Ly, 1.0_dp)/(ny+1)
   A%dx = dx; A%dy = dy
   end procedure
end submodule
