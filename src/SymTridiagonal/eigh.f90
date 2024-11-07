submodule(specialmatrices_symtridiagonal) symtridiagonal_eigenvalue_decomposition
   use stdlib_linalg, only: eye
   use stdlib_linalg_lapack, only: steqr
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   implicit none(type, external)

   character(*), parameter :: this = "symtridiagonal_eigenvalues"
contains

   ! Request for eigenvector calculation.
   elemental character function eigenvectors_task(required)
      logical(lk), intent(in) :: required
      eigenvectors_task = merge("I", "N", required)
   end function eigenvectors_task

   ! Process STEQR output flags.
   elemental subroutine handle_steqr_info(compz, n, ldz, info, err)
      character, intent(in) :: compz
      integer(ilp), intent(in) :: n, ldz, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         ! Success.
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for compz=", compz)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for n=", n)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for D.")
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for E.")
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for Z.")
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for ldz=", ldz)
      case (-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Invalid value for work.")
      case (1:)
         err = linalg_state_type(this, LINALG_ERROR, "steqr failed to converge.")
      case default
         err = linalg_state_type(this, LINALG_INTERNAL_ERROR, "Unknown error returned by setqr.")
      end select
   end subroutine handle_steqr_info

   module procedure eigh_rdp
   ! Local variables.
   type(linalg_state_type) :: err0
   integer(ilp) :: n, ldz, info
   real(dp), allocatable :: work(:), dv(:), ev(:)
   real(dp), target :: vectors_dummy(1, 1)
   real(dp), pointer :: zmat(:, :)
   character :: task_vectors

   ! Matrix size.
   n = A%n

   ! Allocate arrays.
   dv = A%dv; ev = A%ev

   ! Should eigenvectors be computed?
   task_vectors = eigenvectors_task(present(vectors))
   if (present(vectors)) then
      vectors = eye(n); zmat => vectors
      allocate (work(2*n - 2)); ldz = n
   else
      zmat => vectors_dummy; allocate (work(1)); ldz = 1
   end if

   ! Compute eigenvalues and eigenvectors.
   call steqr(task_vectors, n, dv, ev, zmat, ldz, work, info)
   call handle_steqr_info(task_vectors, n, n, info, err0)

   ! Return results.
   lambda = dv
   end procedure

   module procedure eigvalsh_rdp
   call eigh(A, lambda)
   end procedure

end submodule symtridiagonal_eigenvalue_decomposition
