submodule(specialmatrices_bidiagonal) bidiagonal_matvecs
   use stdlib_linalg_lapack, only: lagtm
   implicit none(type, external)
contains

   module procedure spmv
   ! Local variables.
   character :: trans
   integer(ilp) :: n, nrhs, ldx, ldy
   real(dp), parameter :: alpha = 1.0_dp, beta = 0.0_dp
   real(dp), allocatable :: dummy(:)
   real(dp), pointer :: xmat(:, :), ymat(:, :)

   ! Setup variables.
   n = A%n; nrhs = 1; ldx = n; ldy = n; trans = "N"; allocate (dummy(n-1)); dummy = 0.0_dp
   y = x; xmat(1:n, 1:nrhs) => x; ymat(1:n, 1:nrhs) => y
   ! Matrix-vector product.
   select case(A%which)
   case("L")
      call lagtm(trans, n, nrhs, alpha, A%ev, A%dv, dummy, xmat, ldx, beta, ymat, ldy)
   case("U")
      call lagtm(trans, n, nrhs, alpha, dummy, A%dv, A%ev, xmat, ldx, beta, ymat, ldy)
   end select

   end procedure

   module procedure spmvs
   ! Local variables.
   character :: trans
   integer(ilp) :: n, nrhs, ldx, ldy
   real(dp), parameter :: alpha = 1.0_dp, beta = 0.0_dp
   real(dp), target, allocatable :: dummy(:)

   ! Setup variables.
   n = A%n; nrhs = size(x, 2); ldx = n; ldy = n; trans = "N"
   allocate (dummy(n-1)); dummy = 0.0_dp; y = x
   ! Matrix-vector product.
   select case(A%which)
   case("L")
      call lagtm(trans, n, nrhs, alpha, A%ev, A%dv, dummy, x, ldx, beta, y, ldy)
   case("U")
      call lagtm(trans, n, nrhs, alpha, dummy, A%dv, A%ev, x, ldx, beta, y, ldy)
   end select

   end procedure
end submodule
