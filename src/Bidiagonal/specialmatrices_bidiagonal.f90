module specialmatrices_bidiagonal
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
   public :: eig, eigvals

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   !-----------------------------------------------------------------
   !-----     Base types for Symmetric Bidiagonal matrices     -----
   !-----------------------------------------------------------------

   type, public :: Bidiagonal
      !! Base type used to define a `Bidiagonal` matrix of size `[n, n]`
      !! with diagonals given by rank-1 arrays `dv` (size `n`) and `ev`
      !! (size `n-1`).
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: dv(:), ev(:)
      !! Bidiagonal elements of the matrix.
      character :: which
      !! Whether `A` is lower- or upper-bidiagonal.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Bidiagonal
      !! This interface provides different methods to construct a
      !! `Bidiagonal` matrix. Only the non-zero elements of \( A \) are
      !! stored, i.e.
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       d_1  \\
      !!       e_1  &  d_2   \\
      !!             &  \ddots   &  \ddots   \\
      !!             &           &  e_{n-1} &  d_{n}
      !!    \end{bmatrix}.
      !! \]
      !!
      !! if \(A\) is lower-bidiagonal or
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       d_1   &  e_1 \\
      !!             &  \ddots   &  \ddots   \\
      !!             &           &  d_{n-1}  &  e_{n-1}  \\
      !!             &           &           &  d_n
      !!    \end{bmatrix}
      !! \]
      !!
      !! if \(A\) is upper-bidiagonal.
      !!
      !! @warning
      !! By default, the matrix is lower-bidiagonal. To create an upper-
      !! bidiagonal, set `A%which = "U"`.
      !! @endwarning
      !!
      !! #### Syntax
      !!
      !! - Construct a `Bidiagonal` matrix filled with zeros:
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    type(Bidiagonal) :: A
      !!
      !!    A = Bidiagonal(n)
      !! ```
      !!
      !! - Construct a `Bidiagonal` matrix from rank-1 arrays:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), allocatable :: ev(:), dv(:)
      !!    type(Bidiagonal) :: A
      !!    integer :: i
      !!
      !!    dv = [(i, i=1, n)]; ev = [(2*i, i=1, n)]
      !!    A = Bidiagonal(dv, ev)
      !! ```
      !!
      !! - Construct a `Bidiagonal` matrix with constant diagonals:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), parameter :: d = 1.0_dp, e = 2.0_dp
      !!    type(Bidiagonal) :: A
      !!
      !!    A = Bidiagonal(d, e, n)
      !! ```
      !!
      !! @note
      !! Only `double precision` is currently supported for this matrix type.
      !! @endnote
      pure module function initialize(n) result(A)
         !! Construct a `Bidiagonal` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Bidiagonal) :: A
         !! Symmetric Bidiagonal matrix.
      end function

      pure module function construct(dv, ev, which) result(A)
         !! Construct a `Bidiagonal` matrix from the rank-1 arrays `dv`
         !! and `ev`.
         real(dp), intent(in) :: dv(:), ev(:)
         !! Bidiagonal elements of the matrix.
         character, optional, intent(in) :: which
         !! Whether `A` is lower- or upper-diagonal.
         type(Bidiagonal) :: A
         !! Bidiagonal matrix.
      end function

      pure module function construct_constant(d, e, n, which) result(A)
         !! Construct a `Bidiagonal` matrix with constant diagonal elements.
         real(dp), intent(in) :: d, e
         !! Bidiagonal elements of the matrix.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         character, optional, intent(in) :: which
         !! Whether `A` is lower- or upper-bidiagonal.
         type(Bidiagonal) :: A
         !! Symmetric Bidiagonal matrix.
      end function
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a
      !! `Bidiagonal` matrix, both for matrix-vector and matrix-matrix
      !! products. For a matrix-matrix product \( C = AB \), only the matrix
      !! \( A \) has to be a `Bidiagonal` matrix. Both \( B \) and \( C \)
      !! need to be standard Fortran rank-2 arrays. All the underlying
      !! functions are defined as `pure`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    y = matmul(A, x)
      !! ```
      module function spmv(A, x) result(y)
         !! Compute the matrix-vector product \(y = Ax\) for a `Bidiagonal`
         !! matrix \(A\). Both `x` and `y` are rank-1 arrays with the same
         !! kind as `A`.
         type(Bidiagonal), target, intent(in) :: A
         !! Input matrix.
         real(dp), target, intent(in) :: x(:)
         !! Input vector.
         real(dp), target, allocatable :: y(:)
         !! Output vector.
      end function

      pure module function spmvs(A, X) result(Y)
         !! Compute the matrix-matrix product \(Y = Ax\) for a `Bidiagonal`
         !! matrix \(A\) and a dense matrix \(X\) (rank-2 array). \(Y\) is
         !! also a rank-2 array with the same dimensions as \(X\).
         type(Bidiagonal), intent(in) :: A
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
      !! `Bidiagonal` matrix. It also enables to solve a linear system with
      !! multiple right-hand sides.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    x = solve(A, b)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `b` :  Rank-1 or rank-2 array defining the right-hand side(s).
      !!          It is an `intent(in)` argument.
      !!
      !! - `x` :  Solution of the linear system.
      !!          It has the same type and shape as `b`.
      pure module function solve_single_rhs(A, b) result(x)
         !! Solve the linear system \(Ax=b\) where \(A\) is of type
         !! `Bidiagonal` and `b` a standard rank-1 array. The solution
         !! vector `x` has the same dimension and kind as `b`.
         type(Bidiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable, target :: x(:)
         !! Solution vector.
      end function

      pure module function solve_multi_rhs(A, b) result(x)
         !! Solve the linear system \(AX=B\) where \(A\) is of type
         !! `Bidiagonal` and `B` a standard rank-2 array. The solution matrix
         !! `X` has the same dimensions and kind as `B`.
         type(Bidiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable, target :: x(:, :)
         !! Solution vectors.
      end function
   end interface

   interface inv
      pure module function inv_rdp(A) result(B)
         !! Utility function to compute the inverse of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
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
      !! `Bidiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    d = det(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is in an `intent(in)` argument.
      !!
      !! - `d` :  Determinant of the matrix.
      pure module function det_rdp(A) result(d)
         !! Compute the determinant of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: d
         !! Determinant of the matrix.
      end function
   end interface

   interface trace
      !! This interface overloads the `trace` interface from `stdlib_linalg`
      !! to compute the trace of a matrix \( A \) of type `Bidiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    tr = trace(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `tr`:  Trace of the matrix.
      pure module function trace_rdp(A) result(tr)
         !! Compute the trace of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function
   end interface

   !------------------------------------------------
   !-----     Singular Value Decomposition     -----
   !------------------------------------------------

   interface svdvals
      !! This interface overloads the `svdvals` interface from
      !! `stdlib_linalg` to compute the singular values of a `Bidiagonal`
      !! matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    s = svdvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `s` :  Vector of singular values sorted in decreasing order.
      module function svdvals_rdp(A) result(s)
         !! Compute the singular values of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: s(:)
         !! Singular values in descending order.
      end function
   end interface

   interface svd
      !! This interface overloads the `svd` interface from `stdlib_linalg`
      !! to compute the the singular value decomposition of a `Bidiagonal`
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
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `s` :  Rank-1 array `real` array returning the singular values of
      !!          `A`. It is an `intent(out)` argument.
      !!
      !! - `u` (optional) :   Rank-2 array of the same kind as `A` returning
      !!                      the left singular vectors of `A` as columns.
      !!                      Its size should be `[n, n]`. It is an
      !!                      `intent(out)` argument.
      !!
      !! - `vt` (optional):   Rank-2 array of the same kind as `A` returning
      !!                      the right singular vectors of `A` as rows.
      !!                      Its size should be `[n, n]`. It is an
      !!                      `intent(out)` argument.
      module subroutine svd_rdp(A, s, u, vt)
         !! Compute the singular value decomposition of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(out) :: s(:)
         !! Singular values in descending order.
         real(dp), optional, intent(out) :: u(:, :)
         !! Left singular vectors as columns.
         real(dp), optional, intent(out) :: vt(:, :)
         !! Right singular vectors as rows.
      end subroutine
   end interface

   !--------------------------------------------
   !-----     Eigenvalue Decomposition     -----
   !--------------------------------------------

   interface eigvals
      !! This interface overloads the `eigvalsh` interface from
      !! `stdlib_linalg` to compute the eigenvalues of a real-valued matrix
      !! \( A \) whose type is `Bidiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  `real`-valued matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda` :  Vector of eigenvalues in increasing order.
      module function eigvals_rdp(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real
         !! `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         complex(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eig
      !! This interface overloads the `eigh` interface from `stdlib_linalg`
      !! to compute the eigenvalues and eigenvectors of a real-valued matrix
      !! \(A\) whose type is `Bidiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eig(A, lambda [, left] [, right])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  `real`-valued matrix of `Bidiagonal`.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda`  :  Rank-1 `real` array returning the eigenvalues of `A`
      !!                in increasing order. It is an `intent(out)` argument.
      !!
      !! - `left` (optional)  :  `complex` rank-2 array of the same kind as
      !!                         `A` returning the left eigenvectors of `A`.
      !!                         It is an `intent(out)` argument.
      !!
      !! - `right` (optional) :  `complex` rank-2 array of the same kind as
      !!                         `A` returning the right eigenvectors of `A`.
      !!                         It is an `intent(out)` argument.
      !!
      !! @note
      !! No specialized eigensolvers for generic `Bidiagonal` matrices exist
      !! in LAPACK. This routine thus falls back to wrapping the `eig`
      !! procedure from `stdlib_linalg` which uses `*geev` under the hood.
      !! @endnote
      module subroutine eig_rdp(A, lambda, left, right)
         !! Utility function to compute the eigenvalues and eigenvectors of a
         !! `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         complex(dp), intent(out) :: lambda(:)
         !! Eigenvalues.
         complex(dp), optional, intent(out) :: right(:, :), left(:, :)
         !! Eigenvectors.
      end subroutine
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      !! This interface provides methods to convert a `Bidiagonal` matrix
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
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Rank-2 array representation of the matrix \( A \).
      module function dense_rdp(A) result(B)
         !! Utility function to convert a `Bidiagonal` matrix to a
         !! rank-2 array.
         type(Bidiagonal), intent(in) :: A
         !! Input diagonal matrix.
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function
   end interface

   interface transpose
      !! This interface overloads the Fortran `intrinsic` procedure to define
      !! the transpose of a `Bidiagonal` matrix.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = transpose(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Bidiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Resulting transposed matrix. It is of the same type as `A`.
      pure module function transpose_rdp(A) result(B)
         !! Utility function to compute the transpose of a `Bidiagonal`
         !! matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         type(Bidiagonal) :: B
         !! Transpose of the matrix.
      end function
   end interface

   interface size
      pure module function size_rdp(A, dim) result(arr_size)
         !! Utility function to return the size of `Bidiagonal` matrix along
         !! a given dimension.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Size of the matrix along the dimension dim.
      end function
   end interface

   interface shape
      pure module function shape_rdp(A) result(arr_shape)
         !! Utility function to get the shape of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface operator(*)
      pure module function scalar_multiplication_rdp(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a
         !! `Bidiagonal` matrix.
         real(dp), intent(in) :: alpha
         type(Bidiagonal), intent(in) :: A
         type(Bidiagonal) :: B
      end function scalar_multiplication_rdp

      pure module function scalar_multiplication_bis_rdp(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a
         !! `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(Bidiagonal) :: B
      end function scalar_multiplication_bis_rdp
   end interface

end module
