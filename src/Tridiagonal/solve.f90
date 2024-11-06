submodule(specialmatrices_tridiagonal) tridiagonal_linear_solver
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: gtsv, gttrf, gttrs, gtrfs
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   implicit none(type, external)

   character(*), parameter :: this = "tridiagonal_linear_solver"
contains

   ! Process GTTRF
   elemental subroutine handle_gttrf_info(err, info)
      !> Error handler.
      type(linalg_state_type), intent(inout) :: err
      ! GTTRF return flag.
      integer(ilp), intent(in) :: info

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gttrf")
      end select
   end subroutine handle_gttrf_info

   ! Process GTTRS
   elemental subroutine handle_gttrs_info(err, info)
      !> Error handler.
      type(linalg_state_type), intent(inout) :: err
      ! GTTRS return flag.
      integer(ilp), intent(in) :: info

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gttrs")
      end select
   end subroutine handle_gttrs_info

   ! Process GTRFS
   elemental subroutine handle_gtrfs_info(err, info)
      !> Error handler.
      type(linalg_state_type), intent(inout) :: err
      ! GTRFS return flag.
      integer(ilp), intent(in) :: info

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gtrfs")
      end select
   end subroutine handle_gtrfs_info

   function tridiagonal_solver(A, b, refine) result(x)
      type(Tridiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: b(:, :)
      !! Right-hand side vectors.
      logical(lk), intent(in) :: refine
      !! Iterative refinement of the solution?
      real(dp), allocatable :: x(:, :)
      !! Solution vectors.

      ! General LAPACK variables.
      integer(ilp) :: n, nrhs, info
      ! LAPACK variables for LU decomposition.
      real(dp), allocatable :: dl(:), d(:), du(:), du2(:)
      integer(ilp), allocatable :: ipiv(:)
      ! LAPACK variables for iterative refinement.
      real(dp), allocatable :: ferr(:), berr(:), work(:)
      integer(ilp), allocatable :: iwork(:)

      ! Error handler.
      type(linalg_state_type) :: err

      ! Initialize data.
      n = A%n; nrhs = size(b, 2); x = b

      !------------------------------------
      !-----     LU factorization     -----
      !------------------------------------

      ! ----- Allocations -----
      allocate (du2(n - 2), ipiv(n))
      dl = A%dl; d = A%dv; du = A%du; 
      ! ----- LU factorization -----
      call gttrf(n, dl, d, du, du2, ipiv, info)
      call handle_gttrf_info(err, info)

      !-------------------------------------
      !-----     Tridiagonal solve     -----
      !-------------------------------------

      ! ----- Solve the system -----
      call gttrs("N", n, nrhs, dl, d, du, du2, ipiv, x, n, info)
      call handle_gttrs_info(err, info)

      !----------------------------------------
      !-----     Iterative refinement     -----
      !----------------------------------------

      if (refine) then
         ! ----- Allocate arrays -----
         allocate (ferr(nrhs), berr(nrhs), work(3*n), iwork(n))
         ! ----- Refinement step -----
         call gtrfs("N", n, nrhs, A%dl, A%dv, A%du, dl, d, du, du2, ipiv, b, &
                    n, x, n, ferr, berr, work, iwork, info)
      end if
   end function tridiagonal_solver

   module procedure solve_single_rhs
   ! Local variables.
   logical(lk) :: refine_
   real(dp), pointer :: xmat(:, :), bmat(:, :)
   refine_ = optval(refine, .false.)
   x = b; xmat(1:A%n, 1:1) => x; bmat(1:A%n, 1:1) => b
   xmat = tridiagonal_solver(A, bmat, refine_)
   end procedure

   module procedure solve_multi_rhs
   ! Local variables.
   logical(lk) :: refine_
   refine_ = optval(refine, .false.)
   x = tridiagonal_solver(A, b, refine_)
   end procedure
end submodule
