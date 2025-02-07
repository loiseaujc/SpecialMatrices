submodule(specialmatrices_poisson2D) poisson2D_matvecs
   implicit none(type, external)
contains
   module procedure spmv
   real(dp), pointer, contiguous::  xmat(:, :), ymat(:, :)
   real(dp) :: dx, dy, d2x, d2y
   integer(ilp) :: i, j, nx, ny

   ! Initialize arrays and pointers.
   nx = A%nx; ny = A%ny; dx = A%dx; dy = A%dy
   y = x; ymat(1:nx, 1:ny) => y; xmat(1:nx, 1:ny) => x

   ! Evaluate the Laplace operator.
   do concurrent (i=1:nx, j=1:ny)
      ! Horizontl contribution.
      if (i == 1) then
         d2x = (-2*xmat(i, j) + xmat(i+1, j))/dx**2
      else if (i == nx) then
         d2x = (xmat(i-1, j) - 2*xmat(i, j))/dx**2
      else
         d2x = (xmat(i-1, j) -2*xmat(i, j) + xmat(i+1, j))/dx**2
      endif
   
      ! Vertical contribution.
      if (j == 1) then
         d2y = (-2*xmat(i, j) + xmat(i, j+1))/dy**2
      else if (j == ny) then
         d2y = (xmat(i, j-1) - 2*xmat(i, j))/dy**2
      else
         d2y = (xmat(i, j-1) - 2*xmat(i, j) + xmat(i, j+1))/dy**2
      endif

      ! Laplacien.
      ymat(i, j) = d2x + d2y
   enddo 
   end procedure

   module procedure spmvs
   real(dp), pointer, contiguous :: ymat(:, :, :), xmat(:, :, :)
   real(dp) :: dx, dy, d2x, d2y
   integer(ilp) :: i, j, k, nx, ny, nrhs

   ! Initialize arrays and pointers.
   nx = A%nx; ny = A%ny; dx = A%dx; dy = A%dy; nrhs = size(x, 2)
   y = x; ymat(1:nx, 1:ny, 1:nrhs) => y; xmat(1:nx, 1:ny, 1:nrhs) => x

   ! Evaluate the Laplace operator.
   do concurrent (i=1:nx, j=1:ny, k=1:nrhs)
      ! Horizontl contribution.
      if (i == 1) then
         d2x = (-2*xmat(i, j, k) + xmat(i+1, j, k))/dx**2
      else if (i == nx) then
         d2x = (xmat(i-1, j, k) - 2*xmat(i, j, k))/dx**2
      else
         d2x = (xmat(i-1, j, k) -2*xmat(i, j, k) + xmat(i+1, j, k))/dx**2
      endif
   
      ! Vertical contribution.
      if (j == 1) then
         d2y = (-2*xmat(i, j, k) + xmat(i, j+1, k))/dy**2
      else if (j == ny) then
         d2y = (xmat(i, j-1, k) - 2*xmat(i, j, k))/dy**2
      else
         d2y = (xmat(i, j-1, k) - 2*xmat(i, j, k) + xmat(i, j+1, k))/dy**2
      endif

      ! Laplacien.
      ymat(i, j, k) = d2x + d2y
   enddo
   end procedure
end submodule
