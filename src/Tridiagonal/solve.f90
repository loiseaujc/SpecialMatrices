submodule(specialmatrices_tridiagonal) tridiagonal_linear_solver
   use stdlib_linalg_lapack, only: gtsv
   implicit none(type, external)
contains
   module procedure solve_single_rhs
   integer(ilp) :: i, n, nrhs, info
   real(dp), allocatable :: dl(:), d(:), du(:)
   real(dp), pointer :: xmat(:, :)
   ! Initialize arrays.
   n = A%n; nrhs = 1; dl = A%dl; d = A%dv; du = A%du
   x = b; xmat(1:n, 1:nrhs) => x
   call gtsv(n, nrhs, dl, d, du, xmat, n, info)
   end procedure

   module procedure solve_multi_rhs
   integer(ilp) :: i, n, nrhs, info
   real(dp), allocatable :: dl(:), d(:), du(:)
   real(dp), pointer :: xmat(:, :)
   ! Initialize arrays.
   n = A%n; nrhs = size(b, 2); dl = A%dl; d = A%dv; du = A%du
   x = b; xmat(1:n, 1:nrhs) => x
   call gtsv(n, nrhs, dl, d, du, xmat, n, info)
   end procedure
end submodule
