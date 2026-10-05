submodule(specialmatrices_toeplitz) toeplitz_linear_solver
   use stdlib_linalg, only: norm
   use stdlib_linalg_iterative_solvers, only: stdlib_linop_t => stdlib_linop_dp_type, &
                                              stdlib_gmres => stdlib_solve_gmres_kernel, &
                                              workspace_t => stdlib_solver_workspace_dp_type, &
                                              stdlib_size_wksp_gmres

   implicit none(type, external)
contains
   module procedure solve_single_rhs
   integer(ilp) :: m, n
   type(stdlib_linop_t) :: linop, precond
   type(Circulant) :: P

   !> Sanity checks.
   m = size(A, 1); n = size(A, 2)
   if (m /= n) error stop "Matrix is not square."
   if (size(b) /= n) error stop "Dimension of b is inconsistent with A."

   !> Allocate solution vector.
   allocate (x(n), source=0.0_dp)

   !> stdlib_linop wrapper for linop and precond.
   linop%matvec => matvec
   P = strang_preconditioner(A)
   precond%matvec => precond_matvec

   !> Solve the linear system with preconditioned GMRES.
   block
      real(dp), parameter :: atol = epsilon(1.0_dp)
      real(dp), parameter :: rtol = sqrt(atol)
      integer(ilp), parameter :: maxiter = 1024
      integer(ilp), parameter :: kdim = 32
      type(workspace_t) :: workspace
      logical(lk), parameter :: compact = .true.
      allocate (workspace%tmp(n, stdlib_size_wksp_gmres(kdim, compact)), source=0.0_dp)
      call stdlib_gmres(linop, precond, b, x, rtol, atol, maxiter, kdim, workspace, compact)
   end block
contains
   pure subroutine matvec(x, y, alpha, beta, op)
      real(dp), intent(in) :: x(:)
      real(dp), intent(inout) :: y(:)
      real(dp), intent(in) :: alpha, beta
      character(1), intent(in) :: op
      y = alpha*matmul(A, x) + beta*y
   end subroutine matvec

   pure subroutine precond_matvec(x, y, alpha, beta, op)
      real(dp), intent(in) :: x(:), alpha, beta
      real(dp), intent(inout) :: y(:)
      character(1), intent(in) :: op
      y = alpha*matmul(P, x) + beta*y
   end subroutine precond_matvec
   end procedure solve_single_rhs

   module procedure solve_multi_rhs
   integer(ilp) :: i
   allocate (x, mold=b)
   do i = 1, size(b, 2)
      x(:, i) = solve(A, b(:, i))
   end do
   end procedure solve_multi_rhs

   !--------------------------------------------
   !-----     CIRCULANT PRECONDITIONER     -----
   !--------------------------------------------

   pure function strang_preconditioner(T) result(C)
      type(Toeplitz), intent(in) :: T
      type(Circulant)            :: C
      real(dp), allocatable      :: c_vec(:)
      integer(ilp)               :: i, n, n2

      !> Dimension of the matrix.
      n = size(T, 1); n2 = n/2

      !> Circulant vector.
      allocate (c_vec(n))
      do concurrent(i=1:n2 + 1)
         c_vec(i) = T%vc(i)
      end do
      do concurrent(i=n2 + 2:n)
         c_vec(i) = T%vr(n - i + 2)
      end do

      !> Circulant matrix.
      C = Circulant(c_vec)
   end function strang_preconditioner

end submodule toeplitz_linear_solver
