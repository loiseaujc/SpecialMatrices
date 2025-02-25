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
      pure module function initialize(nx, ny, Lx, Ly) result(A)
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
      module function spmv(A, x) result(y)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), target, intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable, target :: y(:)
         !! Output vector.
      end function

      module function spmvs(A, x) result(y)
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
      module function solve_single_rhs(A, b) result(x)
         type(Poisson2D), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in), target :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable, target :: x(:)
         !! Solution vector.
      end function

      module function solve_multi_rhs(A, b) result(x)
         type(Poisson2D), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in), target :: b(:, :)
         !! Right-hand side vector.
         real(dp), allocatable, target :: x(:, :)
         !! Solution vector.
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
      module function eigvalsh_rdp(A) result(lambda)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, target :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eigh
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
      pure module function dense_rdp(A) result(B)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: B(:, :)
         !! Dense representation.
      end function
   end interface

   interface shape
      pure module function shape_rdp(A) result(arr_shape)
         type(Poisson2D), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface size
      pure module function size_rdp(A, dim) result(arr_size)
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
