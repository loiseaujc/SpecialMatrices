module specialmatrices_strang
   use stdlib_linalg_constants, only: dp, ilp, lk
   implicit none(type, external)
   private

   ! --> Linear algebra
   public :: det, trace
   public :: matmul
   public :: solve
   public :: eigh, eigvalsh

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size

   !---------------------------------------------------
   !-----     Base type for the Strang matrix     -----
   !---------------------------------------------------

   type, public :: Strang
      !! Base type used to define a `Strang` matrix.
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Strang
      !! Constructor for generating the `Strang` matrix of size `n`. Only
      !! `double precision` is currently supported. The matrix corresponds
      !! to the standard 3-point finite-difference approximation of the
      !! 1D Laplace operator with unit grid-spacing (\(\Delta x = 1\)) and
      !! homogeneous Dirichlet boundary conditions.
      pure module function initialize(n) result(A)
         integer(ilp), intent(in) :: n
         type(Strang) :: A
      end function
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for the `Strang`
      !! matrix. It is overloaded both for matrix-vector and matrix-matrix products.!! For matrix-matrix product \( C = A B \), only \(A\) can be a `Strang`
      !! matrix. Both \( B \) and \( C \) are standard rank-2 arrays.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    y = matmul(A, x)
      !!    C = matmul(A, B)
      !! ```
      pure module function spmv(A, x) result(y)
         !! Driver for the matrix-vector product.
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function

      pure module function spmvs(A, X) result(Y)
         !! Driver for the matrix-matrix product.
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: X(:, :)
         !! Input vector.
         real(dp), allocatable :: Y(:, :)
         !! Output vector.
      end function
   end interface

   !-----------------------------------------------
   !-----     Linear systems of equations     -----
   !-----------------------------------------------

   interface solve
      !! This interface overloads the `solve` interface from `stdlib_linalg` for solving
      !! a system \(Ax = b\) where \(A\) is a `Strang` matrix. It also enables to
      !! solve a linear system with multiple right-hand sides.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    x = solve(A, b)
      !! ```
      module function solve_single_rhs(A, b, refine) result(x)
         type(Strang), intent(in) :: A
         !! Coefficient matrix.
         real(dp), target, intent(in) :: b(:)
         !! Right-hand side vector.
         logical(lk), optional, intent(in) :: refine
         !! Whether iterative refinement is used or not.
         real(dp), allocatable, target :: x(:)
         !! Solution vector.
      end function

      module function solve_multi_rhs(A, b, refine) result(x)
         type(Strang), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:, :)
         !! Right-hand side vectors.
         logical(lk), optional, intent(in) :: refine
         !! Whether iterative refinement is used or not.
         real(dp), allocatable :: x(:, :)
         !! Solution vectors.
      end function
   end interface

   !-----------------------------------------
   !-----     Determinant and Trace     -----
   !-----------------------------------------

   interface det
      pure module function det_rdp(A) result(d)
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp) :: d
         !! Determinant of the matrix.
      end function
   end interface

   interface trace
      pure module function trace_rdp(A) result(tr)
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function
   end interface

   !--------------------------------------------
   !-----     Eigenvalue Decomposition     -----
   !--------------------------------------------

   interface eigvalsh
      !! This interface overloads the `eigvalsh` interface from `stdlib_linalg`
      !! to compute the eigenvalues of a `Strang` matrix. Note that these
      !! eigenvalues are known analytically.
      !!
      !! ### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvalsh(A)
      !! ```
      pure module function eigvalsh_rdp(A) result(lambda)
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eigh
      !! This interface overloads the `eigh` interface from `stdlib_linalg`
      !! to compute the eigenvalues and eigenvectors of a `Strang` matrix.
      !! Note that both quantities are known analytically.
      !!
      !! ```fortran
      !!    call eigh(A, lambda, vectors)
      !! ```
      pure module subroutine eigh_rdp(A, lambda, vectors)
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out) :: lambda(:)
         !! Eigenvalues.
         real(dp), allocatable, intent(out) :: vectors(:, :)
         !! Eigenvectors.
      end subroutine
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      !! This interface provides a method to convert a matrix of type `Strang`
      !! to its dense representation as a standard rank-2 array
   .  pure module function dense_rdp(A) result(B)
         type(Strang), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: B(:, :)
         !! Dense representation.
      end function
   end interface

   interface shape
      pure module function shape_rdp(A) result(arr_shape)
         type(Strang), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface size
      pure module function size_rdp(A, dim) result(arr_size)
         type(Strang), intent(in) :: A
         !! Input matrix.
         integer(ilp), intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Corresponding size.
      end function
   end interface
contains
end module
