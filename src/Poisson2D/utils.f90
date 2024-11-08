submodule(specialmatrices_poisson2D) poisson2D_utils
   use stdlib_linalg, only: eye, kron => kronecker_product
   use specialmatrices_strang
   implicit none(type, external)
contains

   module procedure dense_rdp
   integer(ilp) :: nx, ny, n
   real(dp)     :: dx, dy
   real(dp), allocatable :: Idx(:, :), Idy(:, :), D2x(:, :), D2y(:, :)

   ! Initialize data.
   nx = A%nx; ny = A%ny; dx = A%dx; dy = A%dy

   ! 1D Laplace operator in each direction.
   D2x = -dense(Strang(nx))/dx**2; D2y = -dense(Strang(ny))/dy**2

   ! Corresponding 2D Laplace operator.
   Idx = eye(nx); Idy = eye(ny)
   !B = kron(Idy, D2x) + kron(D2y, Idx)
   B = kron(D2x, Idy) + kron(Idx, D2y)
   end procedure

   module procedure shape_rdp
      arr_shape = A%nx * A%ny
   end procedure

   module procedure size_rdp
      arr_size = A%nx * A%ny
   end procedure
end submodule
