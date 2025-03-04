module specialmatrices_toeplitz
   use stdlib_linalg_constants, only: dp, ilp, lk
   use stdlib_linalg, only: eig, eigvals, svd, svdvals
   use specialmatrices_circulant
   implicit none(type, external)
   private

   ! --> Linear algebra
   public :: transpose
   public :: matmul
   public :: solve
   public :: svd, svdvals
   public :: eig, eigvals

   ! --> Utility functions.
   public :: Circulant
   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   !---------------------------------------------------
   !-----     Base type for Toeplitz matrices     -----
   !---------------------------------------------------

   type, public :: Toeplitz
      !! Base type to define a `Toeplitz` matrix of size [m x n]. The first
      !! column is given by the vector `vc` while the first row is given by
      !! `vr`.
      private
      integer(ilp) :: m, n
      !! Dimensions of the matrix.
      real(dp), allocatable :: vc(:)
      !! First column of the matrix.
      real(dp), allocatable :: vr(:)
      !! First row of the matrix.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Toeplitz
      !! This interface provides methods to construct `Toeplitz` matrices.
      !! Given a vector `vc` specifying the first column of the matrix and a
      !! vector `vr` specifying its first row, the associated `Toeplitz`
      !! matrix is the following \(m \times n\) matrix
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       t_0      &  t_{-1}      &  \cdots   &  t_{-(n-1)}      \\
      !!       t_1      &  t_0      &  \cdots   &  \vdots   \\
      !!       \vdots   &  \ddots   &  \ddots   &  t_{-1}      \\
      !!       t_{m-1}      &  \cdots   &  t_1      &  t_0
      !!    \end{bmatrix}.
      !! \]
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    integer, parameter :: m = 100, n = 200
      !!    real(dp) :: vc(n), vr(n)
      !!    type(Toeplitz) :: A
      !!
      !!    call random_number(vc) ; call random_number(vr)
      !!    A = Toeplitz(vc, vr)
      !! ```
      !!
      !! @warning
      !! The element \( A_{11} \) is read from the first entry of the vector
      !! `vc`. The first entry of `vr` is not referenced.
      !! @endwarning
      !!
      !! @note
      !! Only `double precision` is currently supported for this matrix type.
      !! @endnote
      pure module function construct(vc, vr) result(A)
         !! Construct a `Toeplitz` matrix from the rank-1 arrays `vc` and `vr`.
         real(dp), intent(in) :: vc(:)
         !! First column of the matrix.
         real(dp), intent(in) :: vr(:)
         !! First row of the matrix.
         type(Toeplitz) :: A
         !! Corresponding Toeplitz matrix.
      end function
   end interface

   interface Circulant
      !! Utility function to embed an m x n `Toeplitz` matrix into an
      !! (m+n) x (m+n) `Circulant` matrix.
      pure module function Toeplitz2Circulant(T) result(C)
         type(Toeplitz), intent(in) :: T
         type(Circulant) :: C
      end function
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a
      !! `Toeplitz` matrix, both for matrix-vector and matrix-matrix products.
      !! For a matrix-matrix product \( C = AB \), only the matrix \( A \)
      !! has to be a `Toeplitz` matrix. Both \( B \) and \( C \) need to be
      !! standard Fortran rank-2 arrays. All the underlying functions are
      !! defined as `pure`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    y = matmul(A, x)
      !! ```
      !!
      !! @note
      !! Matrix-vector products for `Toeplitz` matrices can be efficiently
      !! computed by embedding the `Toeplitz` matrix into a `Circulant` matrix
      !! of size `[m+n x m+n]` and using the Fast Fourier Transform provided
      !! by `fftpack`.
      !! @endnote
      pure module function spmv(A, x) result(y)
         !! Compute the matrix-vector product for a `Toeplitz` matrix \(A\).
         !! Both `x` and `y` are rank-1 arrays with the same kind as `A`.
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function

      pure module function spmvs(A, X) result(Y)
         !! Compute the matrix-matrix product for a `Toeplitz` matrix `A`.
         !! Both `X` and `Y` are rank-2 arrays with the same kind as `A`.
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:, :)
         !! Input matrix.
         real(dp), allocatable :: y(:, :)
         !! Output matrix.
      end function
   end interface

   !-----------------------------------------------
   !-----     Linear systems of equations     -----
   !-----------------------------------------------

   interface solve
      !! This interface overloads the `solve` interface from `stdlib_linalg`
      !! for solving a linear system \( Ax = b \) where \( A \) is a `Toeplitz`
      !! matrix. It also enables to solve a linear system with multiple
      !! right-hand sides.
      !!
      !! #### Syntax
      !!
      !! To solve a system with \( A \) being of type `Toeplitz`:
      !!
      !! ```fortran
      !!    x = solve(A, b)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Toeplitz` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `b` :  Rank-1 or rank-2 array defining the right-hand side(s).
      !!          It is an `intent(in)` argument.
      !!
      !! - `x` :  Solution of the linear system.
      !!          It has the same type and shape as `b`.
      !!
      !! @note
      !! Under the hood, a `gmres` solver is being used along with a
      !! Circulant preconditioner. By design, `gmres` is run until a
      !! relative tolerance of \(10^{-8}\) is reached.
      !! @endnote
      pure module function solve_single_rhs(A, b) result(x)
         !! Solve the linear system \(Ax=b\) where \(A\) is `Toeplitz` and `b`
         !! a standard rank-1 array. The solution vector `x` has the same
         !! dimension and kind as the right-hand side vector `b`.
         type(Toeplitz), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function

      pure module function solve_multi_rhs(A, B) result(X)
         !! Solve the linear system \(AX=B\), where `A` is `Toeplitz` and `B`
         !! is a rank-2 array. The solution matrix `X` has the same dimension
         !! and kind as the right-hand side matrix `B`.
         type(Toeplitz), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function
   end interface

   !------------------------------------------------
   !-----     Singular Value Decomposition     -----
   !------------------------------------------------

   interface svdvals
      !! This interface overloads the `svdvals` interface from `stdlib_linalg`
      !! to compute the singular values of a `Toeplitz` matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    s = svdvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Toeplitz` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `s` :  Vector of singular values sorted in decreasing order.
      !!
      !! @note
      !! No analytic expression exist for the singular values of a general
      !! `Toeplitz` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `svdvals` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module function svdvals_rdp(A) result(s)
         !! Compute the singular values of a `Toeplitz` matrix.
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: s(:)
         !! Singular values in descending order.
      end function
   end interface

   interface svd
      !! This interface overloads the `svd` interface from `stdlib_linalg` to
      !! compute the the singular value decomposition of a `Toeplitz` matrix
      !! \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call svd(A, s, u, vt)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Toeplitz` type.
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
      !!                      size should be `[n, n]`. It is an `intent(out)`
      !!                      argument.
      !!
      !! @note
      !! No analytic expression exist for the singular value of a general
      !! `Toeplitz` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `svdvals` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module subroutine svd_rdp(A, s, u, vt)
         !! Compute the singular value decomposition of a `Toeplitz` matrix.
         type(Toeplitz), intent(in) :: A
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
      !! This interface overloads the `eigvals` interface from `stdlib_linalg`
      !! to compute the eigenvalues of a real-valued matrix \( A \) whose
      !! type is `Toeplitz`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  `real`-valued matrix of `Toeplitz` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda` :  Vector of eigenvalues in increasing order.
      !!
      !! @note
      !! No analytic expression exist for the eigenvalues of a general
      !! `Toeplitz` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `eigvals` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module function eigvals_rdp(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real `Toeplitz`
         !! matrix.
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         complex(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eig
      !! This interface overloads the `eig` interface from `stdlib_linalg` to
      !! compute the eigenvalues and eigenvectors of a real-valued matrix
      !! \(A\) whose type is `Toeplitz`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eig(A, lambda [, left] [, right])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` : `real`-valued matrix of `Toeplitz`.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda`  :  Rank-1 `real` array returning the eigenvalues of `A`
      !!                in increasing order.
      !!                It is an `intent(out)` argument.
      !!
      !! - `left` (optional) :  `complex` rank-2 array of the same kind as `A`
      !!                         returning the left eigenvectors of `A`.
      !!                         It is an `intent(out)` argument.
      !!
      !! - `right` (optional) : `complex` rank-2 array of the same kind as `A`
      !!                         returning the right eigenvectors of `A`.
      !!                         It is an `intent(out)` argument.
      !!
      !! @note
      !! No analytic expression exist for the eigendecomposition of a general
      !! `Toeplitz` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `eig` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module subroutine eig_rdp(A, lambda, left, right)
         !! Utility function to compute the eigenvalues and eigenvectors of a
         !! `Toeplitz` matrix.
         type(Toeplitz), intent(in) :: A
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
      !! Convert a `Toeplitz` matrix to a standard rank-2 array.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = dense(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Toeplitz` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Rank-2 array representation of the matrix \( A \).
      pure module function dense_rdp(A) result(B)
         !! Utility function to convert a `Toeplitz` matrix to a rank-2 array.
         type(Toeplitz), intent(in) :: A
         !! Input diagonal matrix.
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function
   end interface

   interface transpose
      !! This interface overloads the Fortran `intrinsic` procedure to define
      !! the transpose operation of a `Toeplitz` matrix.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = transpose(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Toeplitz` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Resulting transposed matrix. It is of the same type as `A`.
      pure module function transpose_rdp(A) result(B)
         !! Utility function to compute the transpose of a `Toeplitz` matrix.
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         type(Toeplitz) :: B
         !! Transpose of the matrix.
      end function
   end interface

   interface size
      !! Utility function to return the size of `Toeplitz` matrix along a
      !! given dimension.
      pure module function size_rdp(A, dim) result(arr_size)
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Size of the matrix along the dimension dim.
      end function
   end interface

   interface shape
      !! Utility function to return the size of a `Toeplitz` matrix.
      pure module function shape_rdp(A) result(arr_shape)
         !! Utility function to get the shape of a `Toeplitz` matrix.
         type(Toeplitz), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface operator(*)
      pure module function scalar_multiplication_rdp(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a `Toeplitz` matrix.
         real(dp), intent(in) :: alpha
         type(Toeplitz), intent(in) :: A
         type(Toeplitz) :: B
      end function scalar_multiplication_rdp

      pure module function scalar_multiplication_bis_rdp(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a `Toeplitz` matrix.
         type(Toeplitz), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(Toeplitz) :: B
      end function scalar_multiplication_bis_rdp
   end interface
end module
