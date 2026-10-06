module specialmatrices_hankel
   use stdlib_linalg_constants, only: dp, ilp, lk
   use stdlib_linalg, only: eig, eigvals, svd, svdvals
   use specialmatrices_toeplitz, only: Toeplitz, matmul
   implicit none(type, external)
   private

   ! --> Linear algebra
   public :: transpose
   public :: matmul
   public :: solve
   public :: svd, svdvals
   public :: eigh, eigvalsh

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   !---------------------------------------------------
   !-----     Base type for Hankel matrices     -----
   !---------------------------------------------------

   type, public :: Hankel
      !! Base type to define a `Hankel` matrix of size [m x n] generate from
      !! the vector v.
      private
      integer(ilp) :: m, n
      !! Dimensions of the matrix.
      real(dp), allocatable :: v(:)
      !! Generating vector.
   end type Hankel

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Hankel
      !! This interface provides methods to construct `Hankel` matrices.
      !! Given a vector `vc` specifying the first column of the matrix and a
      !! vector `vr` specifying its last row, the associated `Hankel`
      !! matrix is the following \(m \times n\) matrix
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       h_0      &  h_1      &  \cdots   &  h_{(n-1)}      \\
      !!       h_1      &  h_0      &  \cdots   &  \vdots   \\
      !!       \vdots   &  \ddots   &  \ddots   &  h_{n-1+m-2}      \\
      !!       h_{m-1}      &  \cdots   &  t_1      &  h_{n+m-2}
      !!    \end{bmatrix}.
      !! \]
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    integer, parameter :: m = 100, n = 200
      !!    real(dp) :: vc(n), vr(n)
      !!    type(Hankel) :: A
      !!
      !!    call random_number(vc) ; call random_number(vr)
      !!    A = hankel(Hc, vr)
      !! ```
      !!
      !! @warning
      !! The element \( A_{m1} \) is read from the last entry of the vector
      !! `vc`. The first entry of `vr` is not referenced.
      !! @endwarning
      !!
      !! @note
      !! Only `double precision` is currently supported for this matrix type.
      !! @endnote
      pure module function construct(v, m, n) result(A)
         implicit none(type, external)
         !! Construct a `Hankel` matrix from the rank-1 arrays `vc` and `vr`.
         real(dp), intent(in) :: v(:)
         !! Generating vector.
         integer(ilp), intent(in) :: m, n
         !! Dimensions of the matrix.
         type(Hankel) :: A
         !! Corresponding hankel matrix.
      end function construct
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a
      !! `Hankel` matrix, both for matrix-vector and matrix-matrix products.
      !! For a matrix-matrix product \( C = AB \), only the matrix \( A \)
      !! has to be a `Hankel` matrix. Both \( B \) and \( C \) need to be
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
      !! Matrix-vector products for `Hankel` matrices can be efficiently
      !! computed by transforming the matrix into a `Toeplitz` one and
      !! embedding the `Toeplitz` matrix into a `Circulant` matrix of size
      !! `[m+n x m+n]` and using the Fast Fourier Transform provided
      !! by `fftpack`.
      !! @endnote
      pure module function spmv(A, x) result(y)
         !! Compute the matrix-vector product for a `Hankel` matrix \(A\).
         !! Both `x` and `y` are rank-1 arrays with the same kind as `A`.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function spmv

      pure module function spmvs(A, X) result(Y)
         !! Compute the matrix-matrix product for a `Hankel` matrix `A`.
         !! Both `X` and `Y` are rank-2 arrays with the same kind as `A`.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:, :)
         !! Input matrix.
         real(dp), allocatable :: y(:, :)
         !! Output matrix.
      end function spmvs
   end interface

   !-----------------------------------------------
   !-----     Linear systems of equations     -----
   !-----------------------------------------------

   interface solve
      !! This interface overloads the `solve` interface from `stdlib_linalg`
      !! for solving a linear system \( Ax = b \) where \( A \) is a `Hankel`
      !! matrix. It also enables to solve a linear system with multiple
      !! right-hand sides.
      !!
      !! #### Syntax
      !!
      !! To solve a system with \( A \) being of type `Hankel`:
      !!
      !! ```fortran
      !!    x = solve(A, b)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Hankel` type.
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
         !! Solve the linear system \(Ax=b\) where \(A\) is `Hankel` and `b`
         !! a standard rank-1 array. The solution vector `x` has the same
         !! dimension and kind as the right-hand side vector `b`.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function solve_single_rhs

      pure module function solve_multi_rhs(A, B) result(X)
         !! Solve the linear system \(AX=B\), where `A` is `Hankel` and `B`
         !! is a rank-2 array. The solution matrix `X` has the same dimension
         !! and kind as the right-hand side matrix `B`.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function solve_multi_rhs
   end interface

   !------------------------------------------------
   !-----     Singular Value Decomposition     -----
   !------------------------------------------------

   interface svdvals
      !! This interface overloads the `svdvals` interface from `stdlib_linalg`
      !! to compute the singular values of a `Hankel` matrix \(A\).
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    s = svdvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Hankel` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `s` :  Vector of singular values sorted in decreasing order.
      !!
      !! @note
      !! No analytic expression exist for the singular values of a general
      !! `Hankel` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `svdvals` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module function svdvals_rdp(A) result(s)
         !! Compute the singular values of a `hankel` matrix.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: s(:)
         !! Singular values in descending order.
      end function svdvals_rdp
   end interface

   interface svd
      !! This interface overloads the `svd` interface from `stdlib_linalg` to
      !! compute the the singular value decomposition of a `Hankel` matrix
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
      !! - `A` :  Matrix of `Hankel` type.
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
      !! `Hankel` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `svdvals` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module subroutine svd_rdp(A, s, u, vt)
         !! Compute the singular value decomposition of a `Hankel` matrix.
         implicit none(type, external)
         type(hankel), intent(in) :: A
         !! Input matrix.
         real(dp), intent(out) :: s(:)
         !! Singular values in descending order.
         real(dp), optional, intent(out) :: u(:, :)
         !! Left singular vectors as columns.
         real(dp), optional, intent(out) :: vt(:, :)
         !! Right singular vectors as rows.
      end subroutine svd_rdp
   end interface

   !--------------------------------------------
   !-----     Eigenvalue Decomposition     -----
   !--------------------------------------------

   interface eigvalsh
      !! This interface overloads the `eigvals` interface from `stdlib_linalg`
      !! to compute the eigenvalues of a real-valued matrix \( A \) whose
      !! type is `Hankel`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  `real`-valued matrix of `Hankel` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda` :  Vector of eigenvalues in increasing order.
      !!
      !! @note
      !! No analytic expression exist for the eigenvalues of a general
      !! `Hankel` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `eigvals` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module function eigvalsh_rdp(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real `Hankel`
         !! matrix.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function eigvalsh_rdp
   end interface

   interface eigh
      !! This interface overloads the `eigh` interface from `stdlib_linalg` to
      !! compute the eigenvalues and eigenvectors of a real-valued matrix
      !! \(A\) whose type is `Hankel`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eigh(A, lambda [, left] [, right])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` : `real`-valued matrix of `Hankel`.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda`  :  Rank-1 `real` array returning the eigenvalues of `A`
      !!                in increasing order.
      !!                It is an `intent(out)` argument.
      !!
      !! - `vectors` (optional)  :  `real` rank-2 array of the same kind as `A`
      !!                            returning the left eigenvectors of `A`.
      !!                            It is an `intent(out)` argument.
      !!
      !! @note
      !! No analytic expression exist for the eigendecomposition of a general
      !! `hankel` matrix. Under the hood, the matrix `A` is converted to
      !! its dense representation and the function `eigh` from
      !! `stdlib_linalg` is used.
      !! @endnote
      module subroutine eigh_rdp(A, lambda, vectors)
         !! Utility function to compute the eigenvalues and eigenvectors of a
         !! `Hankel` matrix.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         real(dp), intent(out) :: lambda(:)
         !! Eigenvalues.
         real(dp), optional, intent(out) :: vectors(:, :)
         !! Eigenvectors.
      end subroutine eigh_rdp
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      !! Convert a `Hankel` matrix to a standard rank-2 array.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = dense(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Hankel` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Rank-2 array representation of the matrix \( A \).
      pure module function dense_rdp(A) result(B)
         !! Utility function to convert a `Hankel` matrix to a rank-2 array.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input diagonal matrix.
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function dense_rdp
   end interface

   interface transpose
      !! This interface overloads the Fortran `intrinsic` procedure to define
      !! the transpose operation of a `Hankel` matrix.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = transpose(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Hankel` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Resulting transposed matrix. It is of the same type as `A`.
      pure module function transpose_rdp(A) result(B)
         !! Utility function to compute the transpose of a `hankel` matrix.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         type(Hankel) :: B
         !! Transpose of the matrix.
      end function transpose_rdp
   end interface

   interface size
      !! Utility function to return the size of `Hankel` matrix along a
      !! given dimension.
      pure module function size_rdp(A, dim) result(arr_size)
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Size of the matrix along the dimension dim.
      end function size_rdp
   end interface

   interface shape
      !! Utility function to return the size of a `Hankel` matrix.
      pure module function shape_rdp(A) result(arr_shape)
         !! Utility function to get the shape of a `Hankel` matrix.
         implicit none(type, external)
         type(hankel), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function shape_rdp
   end interface

   interface operator(*)
      pure module function scalar_multiplication_rdp(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a `Hankel` matrix.
         implicit none(type, external)
         real(dp), intent(in) :: alpha
         type(Hankel), intent(in) :: A
         type(Hankel) :: B
      end function scalar_multiplication_rdp

      pure module function scalar_multiplication_bis_rdp(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a `Hankel` matrix.
         implicit none(type, external)
         type(Hankel), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(Hankel) :: B
      end function scalar_multiplication_bis_rdp
   end interface
end module specialmatrices_hankel
