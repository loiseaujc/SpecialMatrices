submodule(specialmatrices_bidiagonal) bidiagonal_linear_solver
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)
contains
   module procedure solve_single_rhs
   integer(ilp) :: i, n, nrhs, info
   real(dp), allocatable :: dl(:), dv(:), du(:)
   real(dp), pointer :: xmat(:, :)
   ! Initialize arrays.
   n = A%n; nrhs = 1
   x = b; xmat(1:n, 1:nrhs) => x
   ! Dispatch to solver.
   select case (A%which)
   case ("L")
      dv = A%dv; dl = A%ev; du = 0.0_dp*A%ev
   case ("U")
      dv = A%dv; dl = 0.0_dp*A%ev; du = A%ev
   end select
   ! Solve.
   call gtsv(n, nrhs, dl, dv, du, xmat, n, info)
   end procedure

   module procedure solve_multi_rhs
   integer(ilp) :: i, n, nrhs, info
   real(dp), allocatable :: dl(:), dv(:), du(:)
   ! Initialize arrays.
   n = A%n; nrhs = size(b, 2); x = b
   ! Dispatch to solver.
   select case (A%which)
   case ("L")
      dv = A%dv; dl = A%ev; du = 0.0_dp*A%ev
   case ("U")
      dv = A%dv; dl = 0.0_dp*A%ev; du = A%ev
   end select
   ! Solve.
   call gtsv(n, nrhs, dl, dv, du, x, n, info)
   end procedure
end submodule
