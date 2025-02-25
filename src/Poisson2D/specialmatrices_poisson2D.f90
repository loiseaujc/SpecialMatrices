module specialmatrices_poisson2D
   use stdlib_linalg_constants, only: dp, ilp, lk
   use specialmatrices_strang, only: Strang, dense, eigvalsh, eigh
   implicit none(type, external)
   private

   public :: Poisson2D
   ! --> Linear algebra.
   public :: det, trace
   public :: matmul
   public :: solve
   public :: eigh, eigvalsh

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size


   !------------------------------------------------------
   !-----     Base type for the Poisson2D matrix     -----
   !------------------------------------------------------

   type :: Poisson2D
      !! Base type used to define a `Poisson2D` matrix on a rectangular domain discretized
      !! with `nx` and `ny` points in each direction and corresponding grid spacings `dx`
      !! and `dy`.
      private
      integer(ilp) :: nx, ny
      !! Dimension of the grid in each direction.
      real(dp) :: dx, dy
      !! Grid spacing in each direction.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Poisson2D
      !! Constructor for generating a `Poisson2D` matrix. Only `double precision` is
      !! supported currently. The matrix corresponds to the standard second-order
      !! accurate finite-difference approximation of the Laplace operator with homogeneous
      !! Dirichlet boundary conditions.
      !!
      !! ### Syntax
      !!
      !! - Construct the finite-difference approximation of the Laplace operator
      !!   on the rectangular domain \(\Omega = \left[ 0, 1 \right] \times \left[ 0, 2 \right]\)
      !!   using 128 points in the horizontal direction and 256 in the vertical one.
      !!
      !! ```fortran
      !!    type(Poisson2D)     :: A
      !!    integer, parameter  :: nx = 128, ny = 256
      !!    real(dp), parameter :: Lx = 1.0_dp, Ly = 2.0_dp
      !!
      !!    A = Poisson2D(nx, ny, Lx, Ly)
      !! ```
      !!
      !! Note that `Lx` and `Ly` are optional. If not specified, they default to `1.0_dp`.
      pure module function initialize(nx, ny, Lx, Ly) result(A)
         !! Utility function to construct a `Poisson2D` matrix.
         integer(ilp), intent(in) :: nx, ny
         !! Number of grid points in each direction.
         real(dp), optional, intent(in) :: Lx, Ly
         !! Physical extent of each dimension.
         type(Poisson2D) :: A
         !! Corresponding Poisson2D matrix.
      end function
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a `Poisson2D` matrix.
      !! The intrinsic `matmul` is overloaded both for matrix-vector and matrix-matrix products.
      !! For a matrix-matrix product \(C = AB\), only the matrix \( A \) has to be a `Poisson2D`
      !! matrix. Both \( B \) and \( C \) are standard Fortran rank-2 arrays.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    y = matmul(A, x)
      !! ```
      module function spmv(A, x) result(y)
         !! Compute the matrix-vector product \(y = Ax\) for a `Poisson2D` matrix \( A \).
         !! Both `x` and `y` are rank-1 arrays with the same kind as `A`.
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), target, intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable, target :: y(:)
         !! Output vector.
      end function

      module function spmvs(A, x) result(y)
         !! Compute the matrix-matrix product \(Y=AX\) for a `Poisson2D` matrix \( A \).
         !! \(X\) and \(Y\) are rank-2 arrays of appropriate size with the same kind as \(A\).
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), target, contiguous, intent(in) :: x(:, :)
         !! Input vectors.
         real(dp), allocatable, target :: y(:, :)
         !! Output vectors.
      end function
   end interface

   !-----------------------------------------------
   !-----     Linear systems of equations     -----
   !-----------------------------------------------

   interface solve
      !! This interface overloads the `solve` interface from `stdlib_linalg` for solving
      !! a system \(Ax = b\) where \(A\) is a `Poisson2D` matrix. It also enables to
      !! solve a linear system with multiple right-hand sides. It uses a fast Poisson solver
      !! leveraging the discrete sine transform provided by `fftpack`. Note that only
      !! homogeneous Dirichlet boundary conditions are handled by default. If non-homogeneous
      !! Dirichlet boundary conditions need to be used, they can be implemented by modifiying
      !! the right-hand side vector. Neuman-type boundary conditions are not supported at all.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    x = solve(A, b)
      !! ```
      pure module function solve_single_rhs(A, b) result(x)
         !! Solve the linear system \(Ax = b\) using a fast Poisson solver.
         type(Poisson2D), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable, target :: x(:)
         !! Solution vector.
      end function

      pure module function solve_multi_rhs(A, b) result(x)
         !! Solve the linear system \(AX=B\) using a fast Poisson solver.
         type(Poisson2D), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable, target :: x(:, :)
         !! Solution vectors.
      end function
   end interface

   !-----------------------------------------
   !-----     Determinant and Trace     -----
   !-----------------------------------------

   interface det
      pure module function det_rdp(A) result(d)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp) :: d
         !! Determinant of the matrix.
      end function
   end interface

   interface trace
      pure module function trace_rdp(A) result(tr)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function
   end interface

   !--------------------------------------------
   !-----     Eigenvalue Decomposition     -----
   !--------------------------------------------

   interface eigvalsh
      !! This interface overloads the `eigvalsh` interface from `stdlib_linalg` to compute the
      !! eigenvalues of the `Poisson2D` matrix \(A\). Note that these eigenvalues are known
      !! analytically.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvalsh(A)
      !! ```
      module function eigvalsh_rdp(A) result(lambda)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, target :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eigh
      !! This interface overloads the `eigh` interface from `stdlib_linalg` to compute the
      !! eigenvalues and eigenvectors of the `Poisson2D` matrix \(A\). Note that both eigenvalues
      !! and eigenvectors are known analytically.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    call eigh(A, lambda, vectors)
      !! ```
      module subroutine eigh_rdp(A, lambda, vectors)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out), target :: lambda(:)
         !! Eigenvalues.
         real(dp), allocatable, intent(out) :: vectors(:, :)
      end subroutine
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      !! This interface provides a method to convert a matrix of type `Poisson2D` to its dense
      !! representation as a standard rank-2 array.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    Adense = dense(A)
      !! ```
      pure module function dense_rdp(A) result(B)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: B(:, :)
         !! Dense representation.
      end function
   end interface

   interface shape
      pure module function shape_rdp(A) result(arr_shape)
         !! Utility function to get the shape of the `Poisson2D` matrix.
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface size
      pure module function size_rdp(A, dim) result(arr_size)
         !! Utility function to return the size of a `Poisson2D` matrix along a given dimension.
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Corresponding size.
      end function
   end interface

contains
end module
