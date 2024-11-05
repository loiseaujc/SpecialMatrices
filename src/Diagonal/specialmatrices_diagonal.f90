module specialmatrices_diagonal
   use stdlib_linalg_constants, only: dp, ilp
   implicit none
   private

   ! --> Linear Algebra.
   public :: transpose
   public :: det, trace
   public :: matmul
   public :: inv
   public :: solve
   public :: svd, svdvals
   public :: eigh, eigvalsh

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   !----------------------------------------------------
   !-----     Base types for Diagonal matrices     -----
   !----------------------------------------------------

   type, public :: Diagonal
      !! Base type used to define a `Diagonal` matrix of size [n x n] with diagonal given by the
      !! rank-1 array `dv`.
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: dv(:)
      !! Diagonal elements of the matrix.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Diagonal
      !! This interface provides different methods to construct a `Diagonal` matrix.
      !! Only `double precision` is supported currently. Only the diagonal elements
      !! of \( A \) are being stored, i.e.
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       d_1 \\
      !!          &  d_2 \\
      !!          &     &  \ddots   \\
      !!          &     &        &  d_n
      !!    \end{bmatrix}.
      !! \]
      !!
      !! ### Syntax
      !!
      !! - Construct a `Diagonal` matrix filled with zeros:
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    type(Diagonal) :: A
      !!
      !!    A = Diagonal(n)
      !! ```
      !!
      !! - Construct a `Diagonal` matrix from a vector.
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    real(dp), allocatable :: dv(:)
      !!    type(Diagonal) :: A
      !!    integer :: i
      !!
      !!    dv = [(i, i=1, n)]; A = Diagonal(dv)
      !! ```
      !!
      !! - Construct a `Diagonal` matrix with constant diagonal element.
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    real(dp), parameter :: d = 2.0_dp
      !!    type(Diagonal) :: A
      !!
      !!    A = Diagonal(d, n)
      !! ```
      !!
      !! - Construct a `Diagonal` matrix from a standard Fortran rank-2 array.
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    real(dp) :: B(n, n)
      !!    type(Diagonal) :: A
      !!
      !!    call random_number(B); A = Diagonal(B)
      !! ```
      pure module function initialize(n) result(A)
         !! Utility function to construct a `Diagonal` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function initialize

      pure module function construct(dv) result(A)
         !! Utility function to construct a `Diagonal` matrix from a rank-1 array.
         real(dp), intent(in) :: dv(:)
         !! Diagonal elements of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function construct

      pure module function construct_constant(d, n) result(A)
         !! Utility function to construct a `Diagonal` matrix with constant diagonal element.
         real(dp), intent(in) :: d
         !! Constant diagonal element of the matrix.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function construct_constant

      module function dense_to_diag(A) result(B)
         !! Utility function to construct a `Diagonal` matrix from a rank-2 array.
         real(dp), intent(in) :: A(:, :)
         !! Dense [n x n] matrix from which to construct the `Diagonal` one.
         type(Diagonal) :: B
         !! Corresponding diagonal matrix.
      end function dense_to_diag
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a `Diagonal` matrix.
      !! The intrinsic `matmul` is overloaded both for matrix-vector and matrix-matrix products.
      !! For a matrix-matrix product \( C = AB \), only the matrix \( A \) has to be a `Diagonal`
      !! matrix. Both \( B \) and \( C \) need to be standard Fortran rank-2 arrays.
      !! All the underlying functions are defined as `pure`.
      !!
      !! #### Syntax
      !!
      !! - For matrix-vector product with `A` being of type `Diagonal` and `x` a standard
      !! rank-1 array:
      !! ```fortran
      !!    y = matmul(A, x)
      !! ```
      !!
      !! - For matrix-matrix product with `A` being of type `Diagonal` and `B` a rank-2 array:
      !! ```fortran
      !!    C = matmul(A, B)
      !! ```
      pure module function spmv(A, x) result(y)
         !! Compute the matrix-vector product \(y = Ax\) for a `Diagonal` matrix \(A\).
         !! Both `x` and `y` are rank-1 arrays with the same kind as `A`.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function

      pure module function spmvs(A, x) result(y)
         !! Compute the matrix-matrix product \(Y = Ax\) for a `Diagonal` matrix \(A\) and a
         !! dense matrix \(X\) (rank-2 array). \(Y\) is also a rank-2 array with the same
         !! dimensions as \(X\).
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: X(:, :)
         !! Input vectors.
         real(dp), allocatable :: Y(:, :)
         !! Output vectors.
      end function
   end interface

   !-----------------------------------------------
   !-----     Linear systems of equations     -----
   !-----------------------------------------------

   interface solve
      !! This interface overloads the `solve` interface from `stdlib_linalg` for solving a linear
      !! system \( Ax = b \) where \( A \) is a `Diagonal` matrix.
      !! It also enables to solve a linear system with multiple right-hand sides.
      !!
      !! #### Syntax
      !!
      !! To solve a system with \( A \) being of type `Diagonal`:
      !!
      !! ```fortran
      !!    x = solve(A, b)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal` type. It is an `intent(in)` argument.
      !!
      !! `b`   :  Rank-1 or rank-2 array defining the right-hand side(s). It is an `intent(in)`
      !!          argument.
      !!
      !! `x`   :  Solution of the linear system. It has the same type and shape as `b`.
      pure module function solve_single_rhs(A, b) result(x)
         !! Solve the linear system \(Ax=b\) where \(A\) is of type `Diagonal` and `b` a
         !! standard rank-1 array. The solution vector `x` has the same dimension and kind
         !! as the right-hand side vector `b`.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function

      pure module function solve_multi_rhs(A, b) result(x)
         !! Solve the linear system \(AX=B\) where \(A\) is of type `Diagonal` and `B` a
         !! standard rank-2 array. The solution matrix `X` has the same dimensions and kind
         !! as the right-hand side matrix `B`.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function
   end interface

   interface inv
      pure module function inv_rdp(A) result(B)
         !! Utility function to compute the inverse of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         type(Diagonal) :: B
         !! Inverse of `A`.
      end function
   end interface

   !-----------------------------------------
   !-----     Determinant and Trace     -----
   !-----------------------------------------

   interface det
      !! This interface overloads the `det` interface from `stdlib_linag` to compute the
      !! determinant \(\det(A)\) where \(A\) is of type `Diagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    d = det(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal` type.
      !!          It is in an `intent(in)` argument.
      !!
      !! `d`   :  Determinant of the matrix.
      pure module function det_rdp(A) result(d)
         !! Compute the determinant of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: d
         !! Determinant of the matrix.
      end function
   end interface

   interface trace
      !! This interface overloads the `trace` interface from `stdlib_linalg` to compute the trace
      !! of a matrix \( A \) of type `Diagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    tr = trace(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `tr`  :  Trace of the matrix.
      pure module function trace_rdp(A) result(tr)
         !! Compute the trace of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function
   end interface

   !------------------------------------------------
   !-----     Singular Value Decomposition     -----
   !------------------------------------------------

   interface svdvals
      !! This interface overloads the `svdvals` interface from `stdlib_linalg` to compute the
      !! singular values of a `Diagonal` matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    s = svdvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `s`   :  Vector of singular values sorted in decreasing order.
      pure module function svdvals_rdp(A) result(s)
         !! Compute the singular values of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: s(:)
         !! Singular values in descending order.
      end function
   end interface

   interface svd
      !! This interface overloads the `svd` interface from `stdlib_linalg` to compute the
      !! the singular value decomposition of a `Diagonal` matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call svd(A, s, u, vt)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `s`   :  Rank-1 array `real` array returning the singular values of `A`.
      !!          It is an `intent(out)` argument.
      !!
      !! `u` (optional) :  Rank-2 array of the same kind as `A` returning the left singular
      !!                   vectors of `A` as columns. Its size should be `[n, n]`.
      !!                   It is an `intent(out)` argument.
      !!
      !! `vt (optional) :  Rank-2 array of the same kind as `A` returning the right singular
      !!                   vectors of `A` as rows. Its size should be `[n, n]`.
      !!                   It is an `intent(out)` argument.
      module subroutine svd_rdp(A, u, s, vt)
         !! Compute the singular value decomposition of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out) :: s(:)
         !! Singular values in descending order.
         real(dp), allocatable, optional, intent(out) :: u(:, :)
         !! Left singular vectors as columns.
         real(dp), allocatable, optional, intent(out) :: vt(:, :)
         !! Right singular vectors as rows.
      end subroutine
   end interface

   !--------------------------------------------
   !-----     Eigenvalue Decomposition     -----
   !--------------------------------------------

   interface eigvalsh
      !! This interface overloads the `eigvalsh` interface from `stdlib_linalg` to compute the
      !! eigenvalues of a real-valued matrix \( A \) whose type is `Diagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvalsh(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  `real`-valued matrix of `Diagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `lambda` :  Vector of eigenvalues in increasing order.
      module function eigvalsh_rdp(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eigh
      !! This interface overloads the `eigh` interface from `stdlib_linalg` to compute the
      !! eigenvalues and eigenvectors of a real-valued matrix \(A\) whose type is `Diagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eigh(A, lambda [, vectors])
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   : `real`-valued matrix of `Diagonal`.
      !!          It is an `intent(in)` argument.
      !!
      !! `lambda` :  Rank-1 `real` array returning the eigenvalues of `A` in increasing order.
      !!             It is an `intent(out)` argument.
      !!
      !! `vectors` (optional) :  Rank-2 array of the same kind as `A` returning the eigenvectors
      !!                         of `A`. It is an `intent(out)` argument.
      module subroutine eigh_rdp(A, lambda, vectors)
         !! Utility function to compute the eigenvalues and eigenvectors of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out) :: lambda(:)
         !! Eigenvalues.
         real(dp), allocatable, optional, intent(out) :: vectors(:, :)
         !! Eigenvectors.
      end subroutine
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      !! This interface provides methods to convert a matrix of one the types defined in
      !! `SpecialMatrix` to a regular rank-2 array.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = dense(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `B`   :  Rank-2 array representation of the matrix \( A \).
      module function dense_rdp(A) result(B)
         !! Utility function to convert a `Diagonal` matrix to a rank-2 array.
         type(Diagonal), intent(in) :: A
         !! Input diagonal matrix.
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function
   end interface

   interface transpose
      !! This interface overloads the Fortran `intrinsic` procedure to define the transpose
      !! operation for the various types defined in `SpecialMatrices`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = transpose(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `B`   :  Resulting transposed matrix. It is of the same type as `A`.
      pure module function transpose_rdp(A) result(B)
         !! Utility function to compute the transpose of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         type(Diagonal) :: B
         !! Transpose of the matrix.
      end function
   end interface

   interface size
      pure module function size_rdp(A, dim) result(arr_size)
         !! Utility function to return the size of `Diagonal` matrix along a given dimension.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Size of the matrix along the dimension dim.
      end function
   end interface

   interface shape
      pure module function shape_rdp(A) result(arr_shape)
         !! Utility function to get the shape of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface operator(*)
      pure module function scalar_multiplication_rdp(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a `Diagonal` matrix.
         real(dp), intent(in) :: alpha
         type(Diagonal), intent(in) :: A
         type(Diagonal) :: B
      end function scalar_multiplication_rdp

      pure module function scalar_multiplication_bis_rdp(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(Diagonal) :: B
      end function scalar_multiplication_bis_rdp
   end interface
end module
