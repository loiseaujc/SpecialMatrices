submodule(specialmatrices_poisson2D) poisson2D_solve
   use fftpack, only: dst => dsint, init_dst => dsinti
   implicit none(type, external)

   character(len=*), parameter :: this = "poisson2D_linear_solver"
contains
   module procedure solve_single_rhs
      ! Local variables.
      integer(ilp) :: nx, ny, i, j, nx_wrk, ny_wrk
      real(dp), allocatable :: lambda_x(:), lambda_y(:)
      real(dp), allocatable :: wsave_x(:), wsave_y(:)
      real(dp), pointer     :: xmat(:, :)
      real(dp) :: scale

      ! Initializes pointer and allocatable arrays.
      nx = A%nx ; ny = A%ny
      x = b ; xmat(1:nx, 1:ny) => x
      nx_wrk = int(2.5*nx + 15, kind=ilp) ; ny_wrk = int(2.5*ny + 15, kind=ilp)
      allocate(wsave_x(nx_wrk)) ; allocate(wsave_y(ny_wrk))

      ! Initializes Discrete Sine Transforms (Type-I).
      call init_dst(nx, wsave_x) ; call init_dst(ny, wsave_y)

      ! Compute the DST-I of the right-hand side.
      do concurrent(j=1:ny)
         call dst(nx, xmat(:, j), wsave_x)
      enddo
      do concurrent(i=1:nx)
         call dst(ny, xmat(i, :), wsave_y)
      enddo

      ! Compute the eigenvalues of the 1D Laplacian operators.
      lambda_x = -eigvalsh(Strang(nx)) / A%dx**2
      lambda_y = -eigvalsh(Strang(ny)) / A%dy**2

      ! Compute the solution in spectral space.
      do concurrent(i=1:nx, j=1:ny)
         xmat(i, j) = xmat(i, j) / (lambda_x(i) + lambda_y(j))
      enddo

      ! Inverse DST-I of the solution.
      scale = 1.0_dp / (2*(nx+1) * 2*(ny+1))
      do concurrent(j=1:ny)
         call dst(nx, xmat(:, j), wsave_x)
      enddo
      do concurrent(i=1:nx)
         call dst(ny, xmat(i, :), wsave_y)
      enddo
      xmat = scale * xmat
  end procedure

   module procedure solve_multi_rhs
      integer(ilp) :: i
      allocate(x, mold=b)
      do concurrent(i=1:size(b, 2))
         x(:, i) = solve(A, b(:, i))
      enddo
   end procedure

end submodule
