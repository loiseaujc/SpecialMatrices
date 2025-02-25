submodule(specialmatrices_poisson2D) poisson2D_eigh
   use stdlib_constants, only: pi => pi_dp
   use stdlib_sorting, only: sort
   implicit none(type, external)
contains
   module procedure eigvalsh_rdp
      integer(ilp) :: i, j
      real(dp), pointer :: lambda_mat(:, :)

      !> Eigenvalues for Dx and Dy.
      real(dp), allocatable :: lambda_x(:), lambda_y(:)
      lambda_x = -eigvalsh(Strang(A%nx)) / A%dx**2
      lambda_y = -eigvalsh(Strang(A%ny)) / A%dy**2

      !> Eigenvalues of the 2D Poisson operator.
      allocate(lambda(A%nx*A%ny)) ; lambda_mat(1:A%nx, 1:A%ny) => lambda
      do concurrent(i=1:A%nx, j=1:A%ny)
         lambda_mat(i, j) = lambda_x(i) + lambda_y(j)
      enddo

      !> Sort eigenvalues.
      call sort(lambda)
   end procedure

   module procedure eigh_rdp
   end procedure
end submodule
