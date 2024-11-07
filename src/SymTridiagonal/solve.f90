submodule(specialmatrices_symtridiagonal) symtridiagonal_linear_solver
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: gttrf, gttrs, gtrfs
   use stdlib_linalg_lapack, only: pttrf, pttrs, ptrfs
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   implicit none(type, external)

   character(*), parameter :: this = "symtridiagonal_linear_solver"
contains

   module procedure solve_single_rhs
   ! Local variables.
   logical(lk) :: refine_
   real(dp), pointer :: xmat(:, :), bmat(:, :)
   refine_ = optval(refine, .false.)
   x = b; xmat(1:A%n, 1:1) => x; bmat(1:A%n, 1:1) => b
   if (A%isposdef) then
      xmat = posdef_symtridiagonal_solver(A, bmat, refine_)
   else
      xmat = symtridiagonal_solver(A, bmat, refine_)
   end if
   end procedure

   module procedure solve_multi_rhs
   ! Local variables.
   logical(lk) :: refine_
   refine_ = optval(refine, .false.)
   if (A%isposdef) then
      x = posdef_symtridiagonal_solver(A, b, refine_)
   else
      x = symtridiagonal_solver(A, b, refine_)
   end if
   end procedure

   !---------------------------------------------------
   !-----     Generic (Sym)Tridiagonal Solver     -----
   !---------------------------------------------------

   ! Process GTTRF
   elemental subroutine handle_gttrf_info(n, info, err)
      integer(ilp), intent(in) :: n, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid problem size n=", n)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid size for dl.")
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid size for d.")
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid size for du.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid size for du2.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid size for ipiv.")
      case (1:)
         err = linalg_state_type(this, LINALG_ERROR, "Singular matrix.")
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gttrf")
      end select
   end subroutine handle_gttrf_info

   ! Process GTTRS
   elemental subroutine handle_gttrs_info(trans, n, nrhs, ldb, info, err)
      character, intent(in) :: trans
      integer(ilp), intent(in) :: n, nrhs, ldb, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for trans", trans)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid problem size n=", n)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid number of rhs nrhs=", nrhs)
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid dimensions for dl.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid dimensions for d.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid dimensions for du2.")
      case (-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid dimensions for ipiv.")
      case (-8)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid dimensions for b.")
      case (-9)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldb=", ldb)
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gttrs")
      end select
   end subroutine handle_gttrs_info

   ! Process GTRFS
   elemental subroutine handle_gtrfs_info(trans, n, nrhs, ldb, ldx, info, err)
      character, intent(in) :: trans
      integer(ilp), intent(in) :: n, nrhs, ldb, ldx, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for trans=", trans)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for n=", n)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for nrhs=", nrhs)
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for dl.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for d.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for du.")
      case (-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for dlf.")
      case (-8)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for df.")
      case (-9)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for duf.")
      case (-10)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for du2.")
      case (-11)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ipiv.")
      case (-12)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for b.")
      case (-13)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldb=", ldb)
      case (-14)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for x.")
      case (-15)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldx=", ldx)
      case (-16)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ferr.")
      case (-17)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for berr.")
      case (-18)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for work.")
      case (-19)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for iwork.")
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gtrfs")
      end select
   end subroutine handle_gtrfs_info

   function symtridiagonal_solver(A, b, refine) result(x)
      type(SymTridiagonal), intent(in) :: A
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
      dl = A%ev; d = A%dv; du = A%ev; 
      ! ----- LU factorization -----
      call gttrf(n, dl, d, du, du2, ipiv, info)
      call handle_gttrf_info(n, info, err)

      !-------------------------------------
      !-----     Tridiagonal solve     -----
      !-------------------------------------

      ! ----- Solve the system -----
      call gttrs("N", n, nrhs, dl, d, du, du2, ipiv, x, n, info)
      call handle_gttrs_info("N", n, nrhs, n, info, err)

      !----------------------------------------
      !-----     Iterative refinement     -----
      !----------------------------------------

      if (refine) then
         ! ----- Allocate arrays -----
         allocate (ferr(nrhs), berr(nrhs), work(3*n), iwork(n))
         ! ----- Refinement step -----
         call gtrfs("N", n, nrhs, A%ev, A%dv, A%ev, dl, d, du, du2, ipiv, b, &
                    n, x, n, ferr, berr, work, iwork, info)
         call handle_gtrfs_info("N", n, nrhs, n, n, info, err)
      end if
   end function symtridiagonal_solver

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
      type(SymTridiagonal), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in) :: b(:, :)
      !! Right-hand side vectors.
      logical(lk), intent(in) :: refine
      !! Iterative refinement of the solution?
      real(dp), allocatable :: x(:, :)
      !! Solution vectors.

      ! General LAPACK variables.
      integer(ilp) :: n, nrhs, info
      ! LAPACK variables for LDL^T decomposition.
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
      ev = A%ev; dv = A%dv
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
         call ptrfs(n, nrhs, A%dv, A%ev, dv, ev, b, n, x, n, ferr, berr, work, info)
         call handle_ptrfs_info(n, nrhs, n, n, info, err)
      end if
   end function posdef_symtridiagonal_solver

end submodule
