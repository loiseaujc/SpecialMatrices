submodule(specialmatrices_strang) strang_linear_solver
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: pttrf, pttrs, ptrfs
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   implicit none(type, external)

   character(*), parameter :: this = "strang_linear_solver"
contains
   module procedure solve_single_rhs
   ! Local variables.
   logical(lk) :: refine_
   real(dp), pointer :: xmat(:, :), bmat(:, :)
   refine_ = optval(refine, .false.)
   x = b; xmat(1:A%n, 1:1) => x; bmat(1:A%n, 1:1) => b
   xmat = posdef_symtridiagonal_solver(A, bmat, refine_)
   end procedure

   module procedure solve_multi_rhs
   ! Local variables.
   logical(lk) :: refine_
   refine_ = optval(refine, .false.)
   x = posdef_symtridiagonal_solver(A, b, refine_)
   end procedure

   !-----------------------------------------------------------
   !-----     Positive-definite SymTridiagonal Solver     -----
   !-----------------------------------------------------------

   ! Process PTTRF
   elemental subroutine handle_pttrf_info(n, info, err)
      integer(ilp), intent(in) :: n, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid matrix dimension n=", n)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for D.")
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for E.")
      case (1:)
         err = linalg_state_type(this, LINALG_ERROR, "Matrix could not be factorized.")
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by pttrf")
      end select
   end subroutine handle_pttrf_info

   ! Process PTTRS
   elemental subroutine handle_pttrs_info(n, nrhs, ldb, info, err)
      integer(ilp), intent(in) :: n, nrhs, ldb, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for n=", n)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for nrhs=", nrhs)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for D.")
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for E.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for B.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldb=", ldb)
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by pttrs")
      end select
   end subroutine handle_pttrs_info

   ! Process PTRFS
   elemental subroutine handle_ptrfs_info(n, nrhs, ldb, ldx, info, err)
      integer(ilp), intent(in) :: n, nrhs, ldb, ldx, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for n=", n)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for nrhs=", nrhs)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for D.")
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for E.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for DF.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for EF.")
      case (-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for B.")
      case (-8)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldb=", ldb)
      case (-9)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for X.")
      case (-10)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldx=", ldx)
      case (-11)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ferr.")
      case (-12)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for berr.")
      case (-13)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for work.")
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by ptrfs")
      end select
   end subroutine handle_ptrfs_info

   function posdef_symtridiagonal_solver(A, b, refine) result(x)
      type(Strang), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: b(:, :)
      !! Right-hand side vectors.
      logical(lk), intent(in) :: refine
      !! Iterative refinement of the solution?
      real(dp), allocatable :: x(:, :)
      !! Solution vectors.

      ! General LAPACK variables.
      integer(ilp) :: i, n, nrhs, info
      ! LAPACK variables for LDL^T decomposition.
      real(dp), allocatable :: dv_mat(:), ev_mat(:)
      real(dp), allocatable :: dv(:), ev(:)
      ! LAPACK variables for iterative refinement.
      real(dp), allocatable :: ferr(:), berr(:), work(:)

      ! Error handler.
      type(linalg_state_type) :: err

      ! Initialize data.
      n = A%n; nrhs = size(b, 2); x = b

      !------------------------------------
      !-----     LU factorization     -----
      !------------------------------------

      ! ----- Allocations -----
      dv = [(2, i=1, n)]; ev = [(-1, i=1, n - 1)]
      dv_mat = dv; ev_mat = ev
      ! ----- LDL^T factorization -----
      call pttrf(n, dv, ev, info)
      call handle_pttrf_info(n, info, err)

      !-------------------------------------
      !-----     Tridiagonal solve     -----
      !-------------------------------------

      ! ----- Solve the system -----
      call pttrs(n, nrhs, dv, ev, x, n, info)
      call handle_pttrs_info(n, nrhs, n, info, err)

      !----------------------------------------
      !-----     Iterative refinement     -----
      !----------------------------------------

      if (refine) then
         ! ----- Allocate arrays -----
         allocate (ferr(nrhs), berr(nrhs), work(2*n))
         ! ----- Refinement step -----
         call ptrfs(n, nrhs, dv_mat, ev_mat, dv, ev, b, n, x, n, ferr, berr, work, info)
         call handle_ptrfs_info(n, nrhs, n, n, info, err)
      end if
   end function posdef_symtridiagonal_solver

end submodule

