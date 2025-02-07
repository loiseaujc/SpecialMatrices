submodule(specialmatrices_bidiagonal) bidiagonal_linear_solver
   use stdlib_linalg_lapack, only: gtsv
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   implicit none(type, external)

   character(*), parameter :: this = "bidiagonal_linear_solver"
contains

   elemental subroutine handle_gtsv_info(n, nrhs, ldb, info, err)
      integer(ilp), intent(in) :: n, nrhs, ldb, info
      type(linalg_state_type), intent(out) :: err

      select case (info)
      case (0)
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for n=", n)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for nrhs=", nrhs)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for DL.")
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for D.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for DU.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for b.")
      case (-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldb=", ldb)
      case (1:)
         err = linalg_state_type(this, LINALG_ERROR, "System is singular.")
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by gtsv.")
      end select
      return
   end subroutine handle_gtsv_info

   module procedure solve_single_rhs
   type(linalg_state_type) :: err0
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
   call handle_gtsv_info(n, nrhs, n, info, err0)
   end procedure

   module procedure solve_multi_rhs
   type(linalg_state_type) :: err0
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
   call handle_gtsv_info(n, nrhs, n, info, err0)
   end procedure
end submodule
