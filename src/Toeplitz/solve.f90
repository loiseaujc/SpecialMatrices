submodule(specialmatrices_toeplitz) toeplitz_linear_solver
   use stdlib_linalg, only: norm
   implicit none(type, external)
contains
   module procedure solve_single_rhs
      integer(ilp) :: m, n

      !> Sanity checks.
      m = size(A, 1) ; n = size(A, 2)
      if (m /= n) error stop "Matrix is not square."
      if (size(b) /= n) error stop "Dimension of b is inconsistent with A."

      allocate(x, mold=b)
      !> Solve the linear system with preconditioned GMRES.
      x = gmres(A, b)
   end procedure

   module procedure solve_multi_rhs
      integer(ilp) :: i
      allocate(x, mold=b)
      ! do concurrent(i=1:size(b, 2))
      do i = 1, size(b, 2)
         x(:, i) = solve(A, b(:, i))
      enddo
   end procedure

   !--------------------------------------------
   !-----     CIRCULANT PRECONDITIONER     -----
   !--------------------------------------------

   pure function strang_preconditioner(T) result(C)
      type(Toeplitz), intent(in) :: T
      type(Circulant)            :: C
      real(dp), allocatable      :: c_vec(:)
      integer(ilp)               :: i, n, n2
      
      !> Dimension of the matrix.
      n = size(T, 1) ; n2 = n/2
      
      !> Circulant vector.
      allocate(c_vec(n))
      do concurrent(i=1:n2+1)
         c_vec(i) = T%vc(i)
      enddo
      do concurrent(i=n2+2:n)
         c_vec(i) = T%vr(n-i+2)
      enddo

      !> Circulant matrix.
      C = Circulant(c_vec)
   end function

   !-------------------------------------------
   !-----     ITERATIVE SOLVER: GMRES     -----
   !-------------------------------------------

   pure function gmres(A, b) result(x)
      type(Toeplitz), intent(in) :: A
      !! Coefficient matrix.
      real(dp), intent(in)       :: b(:)
      !! Right-hand side vector.
      real(dp), allocatable      :: x(:)
      !! Solution vector.

      !> Krylov process.
      integer(ilp), parameter :: kmax = 128   !  Dimension of the Krylov subspace.
      real(dp), allocatable :: H(:, :)       !  Upper Hessenberg matrix.
      real(dp), allocatable :: V(:, :)       !  Krylov basis.
      type(Circulant) :: P                   !  Circulant Preconditioner.

      !> Givens rotations.
      real(dp), allocatable :: c(:), s(:)    !  Cosine and sine components.

      !> Miscellaneous.
      real(dp), allocatable :: e(:), y(:)    !  Right-hand side and solution of the lstsq.
      real(dp), parameter :: atol = 10.0_dp ** (-precision(1.0_dp))  !  Absolute tolerance.
      real(dp), parameter :: rtol = 10.0_dp ** (-8)                  !  Relative tolerance.
      real(dp) :: tol                                                !  Actual tolerance.
      real(dp) :: eps, beta
      logical(lk) :: converged
      integer(ilp) :: i, k

      !----- Allocate and initialize working arrays -----
      allocate(x, mold=b)           ;  x = 0.0_dp  !  Solution vector.
      allocate(H(kmax+1, kmax))     ;  H = 0.0_dp  !  Upper Hessenberg matrix.
      allocate(V(size(b), kmax+1))  ;  V = 0.0_dp  !  Krylov basis.
      allocate(e(kmax+1))           ;  e = 0.0_dp  !  Lstsq right-hand side vector.
      allocate(c(kmax))             ;  c = 0.0_dp  !  Cosine component of the Givens rotations.
      allocate(s(kmax))             ;  s = 0.0_dp  !  Sine component of the Givens rotations.

      !> Preconditioner.
      P = strang_preconditioner(A)
   
      !> Set the tolerance.
      tol = atol + norm(b, 2)*rtol

      !----- Generalized Minimum Residual Method -----
      converged = .false.
      do while (.not. converged)
         !> Initialize data.
         H = 0.0_dp ; V = 0.0_dp ; V(:, 1) = b - matmul(A, x)
         e = 0.0_dp ; e(1) = norm(V(:, 1), 2)
         c = 0.0_dp ; s = 0.0_dp

         if (e(1) > tol) then
            V(:, 1) = V(:, 1) / e(1)
         else
            converged = .true.
            exit
         endif

         !> GMRES iterations.
         gmres_step: do k = 1, kmax
            !> Generate new Krylov vector.
            V(:, k+1) = matmul(A, solve(P, V(:, k)))

            !> Orthogonalization step.
            gram_schmidt: do i = 1, k
               H(i, k) = dot_product(V(:, k+1), V(:, i))
               V(:, k+1) = V(:, k+1) - H(i, k)*V(:, i)
            enddo gram_schmidt

            !> Twice is enough.
            do i = 1, k
               eps = dot_product(V(:, k+1), V(:, i))
               V(:, k+1) = V(:, k+1) - eps*V(:, i)
               H(i, k) = H(i, k) + eps
            enddo

            !> Normalize Krylov vector.
            H(k+1, k) = norm(V(:, k+1), 2)
            if (H(k+1, k) > atol) V(:, k+1) = V(:, k+1) / H(k+1, k)

            !> Apply Givens rotations to the Hessenberg matrix.
            call apply_givens_rotation(H(:k+1, k), c(:k), s(:k))
            !> Update the right-hand side vector accordingly.
            e(k+1) = -s(k)*e(k) ; e(k) = c(k)*e(k)

            !> Check convergence.
            if (abs(e(k+1)) < tol) converged = .true.
            if (converged) exit gmres_step
         enddo gmres_step

         !> Solve the least-squares problem.
         k = min(k, kmax) ; y = solve_triangular(H(:k, :k), e(:k))
         !> Update the solution.
         x = x + solve(P, matmul(V(:, :k), y))
      enddo
   end function

  !------------------------------------
  !-----     GIVENS ROTATIONS     -----
  !------------------------------------

  pure function givens_rotation(x) result(g)
    real(dp), intent(in) :: x(2)
    !! Vector whose second component needs to be eliminated.
    real(dp)             :: g(2)
    !! Entries of the Givens rotation matrix.
    g = x / norm(x, 2)
  end function

  pure subroutine apply_givens_rotation(h, c, s)
    real(dp), intent(inout) :: h(:)
    !! k-th column of the Hessenberg matrix.
    real(dp), intent(inout) :: c(:)
    !! Cosine components of the Givens rotations.
    real(dp), intent(inout) :: s(:)
    !! Sine components of the Givens rotations.

    !----- Internal variables -----
    integer(ilp) :: i, k
    real(dp) :: t, g(2)

    !> Size of the column.
    k = size(h)-1

    !> Apply previous Givens rotations to this new column.
    do i = 1, k-1
      t = c(i)*h(i) + s(i)*h(i+1)
      h(i+1) = -s(i)*h(i) + c(i)*h(i+1)
      h(i) = t
    enddo

    !> Compute the sine and cosine components for the next rotation.
    g = givens_rotation([h(k), h(k+1)]) ; c(k) = g(1) ; s(k) = g(2)

    !> Eliminate H(k+1, k).
    h(k) = c(k)*h(k) + s(k)*h(k+1) ; h(k+1) = 0.0_dp
  end subroutine

  !-------------------------------------------
  !-----     Upper Triangular solver     -----
  !-------------------------------------------

  pure function solve_triangular(A, b) result(x)
    real(dp), intent(in) :: A(:, :)
    !! Matrix to invert.
    real(dp), intent(in) :: b(:)
    !! Right-hand side vector.
    real(dp), allocatable :: x(:)
    !! Solution vector.

    !----- Internal variables ------
    integer(ilp) :: i, j, n

    !> Get problem's dimension.
    n = size(A, 1) ; allocate(x, mold=b) ; x = 0.0_dp

    !> Back-substitution algorithm.
    x(n) = b(n) / A(n, n)
    do i = n-1, 1, -1
      x(i) = b(i) - dot_product(A(i, i+1:), x(i+1:))
      x(i) = x(i) / A(i, i)
    enddo
  end function

end submodule
