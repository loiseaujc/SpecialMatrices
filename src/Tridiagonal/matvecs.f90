submodule(specialmatrices_tridiagonal) tridiagonal_matvecs
   use stdlib_linalg_lapack, only: lagtm
   implicit none(type, external)
contains
   module procedure spmv
   ! Local variables.
   character :: trans
   integer(ilp) :: n, nrhs, ldx, ldy
   real(dp), parameter :: alpha = 1.0_dp, beta = 0.0_dp
   real(dp), pointer :: xmat(:, :), ymat(:, :)

   ! Setup variables.
   n = A%n; nrhs = 1; ldx = n; ldy = n; trans = "N"
   y = x; xmat(1:n, 1:nrhs) => x; ymat(1:n, 1:nrhs) => y
   ! Matrix-vector product.
   call lagtm(trans, n, nrhs, alpha, A%dl, A%dv, A%du, xmat, ldx, beta, ymat, ldy)

   end procedure

   module procedure spmvs
   ! Local variables.
   character :: trans
   integer(ilp) :: n, nrhs, ldx, ldy
   real(dp), parameter :: alpha = 1.0_dp, beta = 0.0_dp

   ! Setup variables.
   n = A%n; nrhs = size(x, 2); ldx = n; ldy = n; trans = "N"; y = x
   ! Matrix-vector product.
   call lagtm(trans, n, nrhs, alpha, A%dl, A%dv, A%du, x, ldx, beta, y, ldy)

   end procedure
end submodule
