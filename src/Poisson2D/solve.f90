submodule(specialmatrices_poisson2D) poisson2D_solve
   use stdlib_linalg_lapack, only: trsyl, hseqr
   use stdlib_linalg, only: eye, diag
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   use stdlib_io_npy, only: save_npy
   implicit none(type, external)

   character(len=*), parameter :: this = "poisson2D_linear_solver"
contains
   module procedure solve_single_rhs
      ! Local variables.
      real(dp), allocatable :: D(:), V(:, :)

      ! Compute eigendecomposition of the Poisson 2D.
      call eigh(A, D, V)

      ! Solution.
      x = matmul(transpose(V), b)
      x = matmul(diag(1.0_dp/D), x)
      x = matmul(V, x)
  end procedure

   module procedure solve_multi_rhs
      ! Local variables.
      real(dp), allocatable :: D(:), V(:, :)

      ! Compute eigendecomposition of the Poisson 2D.
      call eigh(A, D, V)

      ! Solution.
      x = matmul(transpose(V), b)
      x = matmul(diag(1.0_dp/D), x)
      x = matmul(V, x)
   end procedure

   !------------------------------------------------------
   !-----     Hessenberg to Schur canonical form     -----
   !------------------------------------------------------

   subroutine hessenberg_to_schur(H, Z)
      real(dp), intent(inout) :: H(:, :)
      !! On entry : Original Hessenberg matrix.
      !! On exit  : Overwritten by the real-Schur form.
      real(dp), allocatable, intent(out) :: Z(:, :)
      !! Orthonormal real-Schur basis.

      ! LAPACK variables.
      type(linalg_state_type) :: err
      character :: job, compz
      integer(ilp) :: n, ilo, ihi, ldh, ldz, lwork, info
      real(dp), allocatable :: wr(:), wi(:), work(:)
      real(dp) :: work_query(1)

      ! Initialize variables.
      job = "S"; compz = "V"
      n = size(H, 1); ldh = n; 
      Z = eye(n, mold=1.0_dp); ldz = size(Z, 1); ilo = 1; ihi = n
      allocate (wr(n)); allocate (wi(n))

      ! Workspace query.
      call hseqr(job, compz, n, ilo, ihi, H, ldh, wr, wi, Z, ldz, work_query, lwork, info)
      call handle_hseqr_info(job, compz, n, ilo, ihi, ldh, ldz, lwork, info, err)

      ! Compute the real Schur form.
      lwork = work_query(1); allocate (work(lwork))
      call hseqr(job, compz, n, ilo, ihi, H, ldh, wr, wi, Z, ldz, work, lwork, info)
      call handle_hseqr_info(job, compz, n, ilo, ihi, ldh, ldz, lwork, info, err)

      return
   end subroutine

   function solve_sylvester(A, B, C) result(X)
      real(dp), intent(in) :: A(:, :), B(:, :), C(:, :)
      !! Matrices defining the Sylvester equation.
      real(dp), allocatable :: X(:, :)
      !! Solution to the Sylvester equation.

      ! LAPACK variables.
      type(linalg_state_type) :: err
      character :: trana, tranb
      integer(ilp) :: isgn, m, n, lda, ldb, ldc, info
      real(dp) :: scale

      ! Initialize variables.
      trana = "N" ; tranb = "N" ; isgn = 1_ilp
      m = size(A, 1) ; lda = m
      n = size(B, 1) ; ldb = n
      ldc = m

      ! Solve the Sylvester equation.
      X = C ; call trsyl(trana, tranb, isgn, m, n, A, lda, B, ldb, X, ldc, scale, info)
      call handle_info_trsyl(trana, tranb, isgn, m, n, lda, ldb, ldc, info, err)
   end function

   !----------------------------------
   !-----     Error Handling     -----
   !----------------------------------

   elemental subroutine handle_hseqr_info(job, compz, n, ilo, ihi, ldh, ldz, lwork, info, err)
      character, intent(in) :: job, compz
      integer(ilp), intent(in) :: n, ilo, ihi, ldh, ldz, lwork, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         err%state = LINALG_SUCCESS
      case (-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Unknown job =", job)
      case (-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Unknown compz =", compz)
      case (-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal problem size n =", n)
      case (-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for ilo =", ilo)
      case (-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for ihi =", ihi)
      case (-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for H.")
      case (-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for ldh =", ldh)
      case (-8)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for wr.")
      case (-9)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for wi.")
      case (-10)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for Z.")
      case (-11)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for ldz =", ldz)
      case (-12)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for work.")
      case (-13)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for lwork.")
      end select
   end subroutine

   elemental subroutine handle_info_trsyl(trana, tranb, isgn, m, n, lda, ldb, ldc, info, err)
      character, intent(in) :: trana, tranb
      integer(ilp), intent(in) :: isgn, m, n, lda, ldb, ldc, info
      type(linalg_state_type), intent(inout) :: err

      select case (info)
      case (0)
         err%state = LINALG_SUCCESS
      case(-1)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Unknown trana =", trana)
      case(-2)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Unknown tranb =", tranb)
      case(-3)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for isgn =", isgn)
      case(-4)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal problem size m =", m)
      case(-5)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal problem size n =", n)
      case(-6)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for A.")
      case(-7)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for lda =", lda)
      case(-8)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for B.")
      case(-9)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for ldb =", ldb)
      case(-10)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for C.")
      case(-11)
         err = linalg_state_type(this, LINALG_VALUE_ERROR, "Illegal value for ldc = ", ldc)
      end select
   end subroutine


end submodule
