submodule(specialmatrices_poisson2D) poisson2D_eigh
   use stdlib_constants, only: pi => pi_dp
   use stdlib_sorting, only: sort, sort_index
   use stdlib_linalg, only: outer_product
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
      integer(ilp) :: i, j, index(A%nx*A%ny)
      integer(ilp) :: nx, ny, n, counter
      real(dp), pointer :: lambda_mat(:, :)
      real(dp), allocatable :: lambda_x(:), lambda_y(:)
      real(dp), allocatable :: vecs_x(:, :), vecs_y(:, :)

      nx = A%nx ; ny = A%ny ; n = nx*ny

      !> Eigenvalues and eigevectors for Dx and Dy.
      call eigh(Strang(A%nx), lambda_x, vecs_x)
      call eigh(Strang(A%ny), lambda_y, vecs_y)
      !> Scale eigenvalues.
      lambda_x = -lambda_x / A%dx**2
      lambda_y = -lambda_y / A%dy**2

      !> Eigenvalues of the 2D Poisson operator.
      allocate(lambda(n)) ; lambda_mat(1:nx, 1:ny) => lambda
      do concurrent(i=1:nx, j=1:ny)
         lambda_mat(i, j) = lambda_x(i) + lambda_y(j)
      enddo

      !> Eigenvectors of the 2D Poisson operator.
      allocate(vectors(n, n)) ; vectors = 0.0_dp ; counter = 1
      do j = 1, ny
         do i = 1, nx
            vectors(:, counter) = pack(outer_product(vecs_x(:, i), vecs_y(:, j)), .true.)
            counter = counter + 1
         enddo
      enddo

      !> Sort eigenvalues and eigenvectors.
      call sort_index(lambda, index) ; vectors = vectors(:, index)
   end procedure
end submodule
