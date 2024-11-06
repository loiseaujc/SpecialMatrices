submodule(specialmatrices_symtridiagonal) symtridiagonal_linear_solver
   use stdlib_linalg_lapack, only: gtsv, ptsv
   implicit none(type, external)
contains
   module procedure solve_single_rhs
   integer(ilp) :: i, n, nrhs, info
   real(dp), allocatable :: dl(:), d(:), du(:)
   real(dp), pointer :: xmat(:, :)
   ! Initialize arrays.
   n = A%n; nrhs = 1; dl = A%ev; d = A%dv; du = A%ev
   x = b; xmat(1:n, 1:nrhs) => x
   if (A%isposdef) then
      call ptsv(n, nrhs, d, du, xmat, n, info)
   else
      call gtsv(n, nrhs, dl, d, du, xmat, n, info)
   end if
   end procedure

   module procedure solve_multi_rhs
   integer(ilp) :: i, n, nrhs, info
   real(dp), allocatable :: dl(:), d(:), du(:)
   ! Initialize arrays.
   n = A%n; nrhs = size(b, 2); dl = A%ev; d = A%dv; du = A%ev; x = b
   if (A%isposdef) then
      call ptsv(n, nrhs, d, du, x, n, info)
   else
      call gtsv(n, nrhs, dl, d, du, x, n, info)
   end if
   end procedure
end submodule
