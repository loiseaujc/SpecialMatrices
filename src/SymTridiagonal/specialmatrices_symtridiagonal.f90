module specialmatrices_symtridiagonal
   use stdlib_linalg_constants, only: dp, ilp, lk
   implicit none(type, external)
   private

   ! --> Linear algebra
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

   !-----------------------------------------------------------------
   !-----     Base types for Symmetric Tridiagonal matrices     -----
   !-----------------------------------------------------------------

   type, public :: SymTridiagonal
      !! Base type used to define a `SymTridiagonal` matrix of size `[n, n]`
      !! with diagonals given by rank-1 arrays `dv` (size `n`) and `ev`
      !! (size `n-1`).
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: dv(:), ev(:)
      !! SymTridiagonal elements of the matrix.
      logical(lk) :: isposdef
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface SymTridiagonal
      !! This interface provides different methods to construct a
      !! `SymTridiagonal` matrix. Only the non-zero elements of \( A \) are
      !! stored, i.e.
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       d_1   &  e_1  \\
      !!       e_1  &  d_2      &  e_2  \\
      !!             &  \ddots   &  \ddots   &  \ddots   \\
      !!             &           &  e_{n-2} &  d_{n-1}  &  e_{n-1} \\
      !!             &           &           &  e_{n-1} &  d_n
      !!    \end{bmatrix}.
      !! \]
      !!
      !! #### Syntax
      !!
      !! - Construct a `SymTridiagonal` matrix filled with zeros:
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    type(SymTridiagonal) :: A
      !!
      !!    A = SymTridiagonal(n)
      !! ```
      !!
      !! - Construct a `SymTridiagonal` matrix from rank-1 arrays:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), allocatable :: ev(:), dv(:)
      !!    type(SymTridiagonal) :: A
      !!    integer :: i
      !!
      !!    dv = [(i, i=1, n)]; ev = [(2*i, i=1, n)]
      !!    A = Tridiagonal(dv, ev)
      !! ```
      !!
      !! - Construct a `SymTridiagonal` matrix with constant diagonals:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), parameter :: d = 1.0_dp, e = 2.0_dp
      !!    type(SymTridiagonal) :: A
      !!
      !!    A = SymTridiagonal(d, e, n)
      !! ```
      !!
      !! @note
      !! Only `double precision` is currently supported for this matrix type.
      !! @endnote
      !!
      !! @note
      !! If \( A \) is known to be symmetric positive definite, it can be
      !! constructed as `A = SymTridiagonal(dv, ev, ifposdef=.true.)`:w

      !! @endnote
      pure module function initialize(n) result(A)
         !! Construct a `SymTridiagonal` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(SymTridiagonal) :: A
         !! Symmetric Tridiagonal matrix.
      end function

      pure module function construct(dv, ev, isposdef) result(A)
         !! Construct a `SymTridiagonal` matrix from the rank-1 arrays
         !! `dv` and `ev`.
         real(dp), intent(in) :: dv(:), ev(:)
         !! SymTridiagonal elements of the matrix.
         logical(lk), optional, intent(in) :: isposdef
         !! Whether `A` is positive-definite or not.
         type(SymTridiagonal) :: A
         !! Symmetric Tridiagonal matrix.
      end function

      pure module function construct_constant(d, e, n, isposdef) result(A)
         !! Construct a `SymTridiagonal` matrix with constant diagonal
         !! elements.
         real(dp), intent(in) :: d, e
         !! SymTridiagonal elements of the matrix.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         logical(lk), optional, intent(in) :: isposdef
         !! Whether `A` is positive-definite or not.
         type(SymTridiagonal) :: A
         !! Symmetric Tridiagonal matrix.
      end function
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a
      !! `SymTridiagonal` matrix, both for matrix-vector and matrix-matrix
      !! products. For a matrix-matrix product \( C = AB \), only the matrix
      !! \( A \) has to be a `SymTridiagonal` matrix. Both \( B \) and \( C \)
      !! need to be standard Fortran rank-2 arrays. All the underlying
      !! functions are defined as `pure`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    y = matmul(A, x)
      !! ```
      module function spmv(A, x) result(y)
         !! Compute the matrix-vector product \(y = Ax\) for a `SymTridiagonal`
         !! matrix \(A\). Both `x` and `y` are rank-1 arrays with the same
         !! kind as `A`.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), target, intent(in) :: x(:)
         !! Input vector.
         real(dp), target, allocatable :: y(:)
         !! Output vector.
      end function

      pure module function spmvs(A, x) result(y)
         !! Compute the matrix-matrix product \(Y = Ax\) for a `SymTridiagonal`
         !! matrix \(A\) and a dense matrix \(X\) (rank-2 array). \(Y\) is
         !! also a rank-2 array with the same dimensions as \(X\).
         type(SymTridiagonal), intent(in) :: A
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
      !! This interface overloads the `solve` interface from `stdlib_linalg`
      !! for solving a linear system \( Ax = b \) where \( A \) is a
      !! `SymTridiagonal` matrix. It also enables to solve a linear system
      !! with multiple right-hand sides.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    x = solve(A, b [, refine])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `b` :  Rank-1 or rank-2 array defining the right-hand side(s).
      !!          It is an `intent(in)` argument.
      !!
      !! - `refine` (optional) : Logical switch to enable solution refinement.
      !!
      !! - `x` :  Solution of the linear system.
      !!          It has the same type and shape as `b`.
      module function solve_single_rhs(A, b, refine) result(x)
         !! Solve the linear system \(Ax=b\) where \(A\) is of type
         !! `SymTridiagonal` and `b` a standard rank-1 array. The solution
         !! vector `x` has the same dimension and kind as `b`.
         type(SymTridiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), target, intent(in) :: b(:)
         !! Right-hand side vector.
         logical(lk), optional, intent(in) :: refine
         !! Whether iterative refinement of the solution is used or not.
         real(dp), allocatable, target :: x(:)
         !! Solution vector.
      end function

      module function solve_multi_rhs(A, b, refine) result(x)
         !! Solve the linear system \(AX=B\) where \(A\) is of type
         !! `SymTridiagonal` and `B` a standard rank-2 array. The solution
         !! matrix `X` has the same dimensions and kind as `B`.
         type(SymTridiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:, :)
         !! Right-hand side vectors.
         logical(lk), optional, intent(in) :: refine
         !! Whether iterative refinement of the solution is used or not.
         real(dp), allocatable :: x(:, :)
         !! Solution vectors.
      end function
   end interface

   interface inv
      pure module function inv_rdp(A) result(B)
         !! Compute the inverse of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: B(:, :)
         !! Inverse of `A`.
      end function
   end interface

   !-----------------------------------------
   !-----     Determinant and Trace     -----
   !-----------------------------------------

   interface det
      !! This interface overloads the `det` interface from `stdlib_linag` to
      !! compute the determinant \(\det(A)\) where \(A\) is of type
      !! `SymTridiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    d = det(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is in an `intent(in)` argument.
      !!
      !! - `d` :  Determinant of the matrix.
      pure module function det_rdp(A) result(d)
         !! Compute the determinant of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: d
         !! Determinant of the matrix.
      end function
   end interface

   interface trace
      !! This interface overloads the `trace` interface from `stdlib_linalg`
      !! to compute the trace of a matrix \( A \) of type `SymTridiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    tr = trace(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `tr`:  Trace of the matrix.
      pure module function trace_rdp(A) result(tr)
         !! Compute the trace of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function
   end interface

   !------------------------------------------------
   !-----     Singular Value Decomposition     -----
   !------------------------------------------------

   interface svdvals
      !! This interface overloads the `svdvals` interface from `stdlib_linalg`
      !! to compute the singular values of a `SymTridiagonal` matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    s = svdvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `s` :  Vector of singular values sorted in decreasing order.
      module function svdvals_rdp(A) result(s)
         !! Compute the singular values of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: s(:)
         !! Singular values in descending order.
      end function
   end interface

   interface svd
      !! This interface overloads the `svd` interface from `stdlib_linalg` to
      !! compute the the singular value decomposition of a `SymTridiagonal`
      !! matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call svd(A, s [, u] [, vt])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `s` :  Rank-1 array `real` array returning the singular values of
      !!          `A`. It is an `intent(out)` argument.
      !!
      !! - `u` (optional)  :  Rank-2 array of the same kind as `A` returning
      !!                      the left singular vectors of `A` as columns. Its
      !!                      size should be `[n, n]`.
      !!                      It is an `intent(out)` argument.
      !!
      !! - `vt` (optional) :  Rank-2 array of the same kind as `A` returning
      !!                      the right singular vectors of `A` as rows. Its
      !!                      size should be `[n, n]`.
      !!                      It is an `intent(out)` argument.
      module subroutine svd_rdp(A, s, u, vt)
         !! Compute the singular value decomposition of a `SymTridiagonal`
         !! matrix.
         type(SymTridiagonal), intent(in) :: A
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
      !! This interface overloads the `eigvalsh` interface from
      !! `stdlib_linalg` to compute the eigenvalues of a matrix \( A \) whose
      !! type is `SymTridiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvalsh(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  `real`-valued matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda` :  Vector of eigenvalues in increasing order.
      module function eigvalsh_rdp(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real
         !! `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eigh
      !! This interface overloads the `eigh` interface from `stdlib_linalg`
      !! to compute the eigenvalues and eigenvectors of a matrix \(A\) whose
      !! type is `SymTridiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eigh(A, lambda [, vectors])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal`.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda`  :  Rank-1 `real` array returning the eigenvalues of `A`
      !!                in increasing order. It is an `intent(out)` argument.
      !!
      !! - `vectors` (optional)  :  Rank-2 array of the same kind as `A`
      !!                            returning the eigenvectors of `A`.
      !!                            It is an `intent(out)` argument.
      module subroutine eigh_rdp(A, lambda, vectors)
         !! Compute the eigenvalues and eigenvectors of a `SymTridiagonal`
         !! matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out) :: lambda(:)
         !! Eigenvalues.
         real(dp), allocatable, optional, target, intent(out) :: vectors(:, :)
         !! Eigenvectors.
      end subroutine
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      !! This interface provides methods to convert a `SymTridiagonal` matrix
      !! to a regular rank-2 array.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = dense(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Rank-2 array representation of the matrix \( A \).
      module function dense_rdp(A) result(B)
         !! Convert a `SymTridiagonal` matrix to a rank-2 array.
         type(SymTridiagonal), intent(in) :: A
         !! Input diagonal matrix.
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function
   end interface

   interface transpose
      !! This interface overloads the Fortran `intrinsic` procedure to define
      !! the transpose operation for a `SymTridiagonal` matrix.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = transpose(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Resulting transposed matrix. It is of the same type as `A`.
      pure module function transpose_rdp(A) result(B)
         !! Compute the transpose of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         type(SymTridiagonal) :: B
         !! Transpose of the matrix.
      end function
   end interface

   interface size
      pure module function size_rdp(A, dim) result(arr_size)
         !! Return the size of `SymTridiagonal` matrix along a given
         !! dimension.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Size of the matrix along the dimension dim.
      end function
   end interface

   interface shape
      pure module function shape_rdp(A) result(arr_shape)
         !! Return the shape of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface operator(*)
      pure module function scalar_multiplication_rdp(alpha, A) result(B)
         !! Scalar multiplication with a `SymTridiagonal` matrix.
         real(dp), intent(in) :: alpha
         type(SymTridiagonal), intent(in) :: A
         type(SymTridiagonal) :: B
      end function scalar_multiplication_rdp

      pure module function scalar_multiplication_bis_rdp(A, alpha) result(B)
         !! Scalar multiplication with a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(SymTridiagonal) :: B
      end function scalar_multiplication_bis_rdp
   end interface

end module
