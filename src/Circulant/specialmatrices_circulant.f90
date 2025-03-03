module specialmatrices_circulant
   use stdlib_linalg_constants, only: dp, ilp, lk
   use fftpack, only: fft, ifft
   implicit none(type, external)
   private

   ! --> Linear algebra
   public :: transpose
   public :: det, trace
   public :: matmul
   public :: inv
   public :: solve
   ! public :: svd, svdvals
   public :: eig, eigvals

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   !----------------------------------------------------
   !-----     Base type for Circulant matrices     -----
   !----------------------------------------------------

   type, public :: Circulant
      !! Base type to define a `Circulant` matrix of size [n x n] with elements
      !! given by the vector \(c = (c[1], c[2], ..., c[n]).\)
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: c(:)
      !! Generating vector.
      complex(dp), allocatable :: c_hat(:)
      !! Fourier Transform of the generating vector.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Circulant
      !! This interface provides methods to construct `Circulant` matrices.
      !! Given a vector \( \mathbf{c} = (c_1, c_2, \cdots, c_n)\),
      !! the associated `Circulant` matrix is the following `[n x n]` matrix
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       c_1      &  c_2      &  \cdots   &  c_n      \\
      !!       c_n      &  c_1      &  \cdots   &  \vdots   \\
      !!       \vdots   &  \ddots   &  \ddots   &  c_2      \\
      !!       c_2      &  \cdots   &  c_n      &  c_1
      !!    \end{bmatrix}.
      !! \]
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    real(dp) :: c(n)
      !!    type(Circulant) :: A
      !!
      !!    call random_number(c)
      !!    A = Circulant(c)
      !! ```
      !!
      !! @note
      !! Only `double precision` is currently supported for this matrix type.
      !! @endnote
      pure module function initialize(n) result(A)
         !! Construct a `Circulant` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Circulant) :: A
         !! Circulant matrix.
      end function

      pure module function construct(c) result(A)
         !! Construct a `Circulant` matrix from the rank-1 array `c`.
         real(dp), intent(in) :: c(:)
         !! Generating vector.
         type(Circulant) :: A
         !! Corresponding Circulant matrix.
      end function
   end interface

   !-------------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix multiplications     -----
   !-------------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for a
      !! `Circulant` matrix, both for matrix-vector and matrix-matrix products.
      !! For a matrix-matrix product \( C = AB \), only the matrix \( A \)
      !! has to be a `Circulant` matrix. Both \( B \) and \( C \) need to be
      !! standard Fortran rank-2 arrays. All the underlying functions are
      !! defined as `pure`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    y = matmul(A, x)
      !! ```
      pure module function spmv(A, x) result(y)
         !! Compute the matrix-vector product for a `Circulant` matrix \(A\).
         !! Both `x` and `y` are rank-1 arrays with the same kind as `A`.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function

      pure module function spmvs(A, X) result(Y)
         !! Compute the matrix-matrix product for a `Circulant` matrix `A`.
         !! Both `X` and `Y` are rank-2 arrays with the same kind as `A`.
         type(Circulant), intent(in) :: A
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
      !! for solving a linear system \( Ax = b \) where \( A \) is a
      !! `Circulant` matrix. It also enables to solve a linear system with
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
      !! - `A` :  Matrix of `Circulant` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `b` :  Rank-1 or rank-2 array defining the right-hand side(s).
      !!          It is an `intent(in)` argument.
      !!
      !! - `x` :  Solution of the linear system.
      !!          It has the same type and shape as `b`.
      !!
      !! @note
      !! Linear systems characterized by a circulant matrix can be solved
      !! efficiently in \(\mathcal{O}(n \log\ n)\) operations using the
      !! Fast Fourier Transform algorithm available via `fftpack`.
      !! @endnote
      pure module function solve_single_rhs(A, b) result(x)
         !! Solve the linear system \(Ax=b\) where \(A\) is `Circulant` and
         !! `b` a standard rank-1 array. The solution vector `x` has the same
         !! dimension and kind as the right-hand side vector `b`.
         type(Circulant), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function

      pure module function solve_multi_rhs(A, B) result(X)
         !! Solve the linear system \(AX=B\), where `A` is `Circulant` and
         !! `B` is a rank-2 array. The solution matrix `X` has the same
         !! dimension and kind as the right-hand side matrix `B`.
         type(Circulant), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function
   end interface

   interface inv
      pure module function inv_rdp(A) result(B)
         !! Utility function to compute the inverse of a `Circulant` matrix.
         !! If `A` is circulant, its inverse also is circulant.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         type(Circulant) :: B
         !! Inverse of `A`.
      end function
   end interface

   !------------------------------------------
   !-----     Determinant and Trace      -----
   !------------------------------------------

   interface det
      !! This interface overloads the `det` interface from `stdlib_linag` to
      !! compute the determinant \(\det(A)\) where \(A\) is of type `Circulant`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    d = det(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Circulant` type.
      !!          It is in an `intent(in)` argument.
      !!
      !! - `d` :  Determinant of the matrix.
      pure module function det_rdp(A) result(d)
         !! Compute the determinant of a `Circulant` matrix.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         real(dp) :: d
         !! Determinant of the matrix.
      end function
   end interface

   interface trace
      !! This interface overloads the `trace` interface from `stdlib_linalg`
      !! to compute the trace of a matrix \( A \) of type `Circulant`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    tr = trace(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Circulant` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `tr`:  Trace of the matrix.
      pure module function trace_rdp(A) result(tr)
         !! Compute the trace of a `Circulant` matrix.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function
   end interface

   !------------------------------------------------
   !-----     Singular Value Decomposition     -----
   !------------------------------------------------

   ! interface svdvals
   !    !! This interface overloads the `svdvals` interface from `stdlib_linalg` to compute the
   !    !! singular values of a `Circulant` matrix \(A\).
   !    !!
   !    !! #### Syntax
   !    !!
   !    !! ```fortran
   !    !!    s = svdvals(A)
   !    !! ```
   !    !!
   !    !! #### Arguments
   !    !!
   !    !! `A`   :  Matrix of `Circulant` type.
   !    !!          It is an `intent(in)` argument.
   !    !!
   !    !! `s`   :  Vector of singular values sorted in decreasing order.
   !    module function svdvals_rdp(A) result(s)
   !       !! Compute the singular values of a `Circulant` matrix.
   !       type(Circulant), intent(in) :: A
   !       !! Input matrix.
   !       real(dp), allocatable :: s(:)
   !       !! Singular values in descending order.
   !    end function
   ! end interface
   !
   ! interface svd
   !    !! This interface overloads the `svd` interface from `stdlib_linalg` to compute the
   !    !! the singular value decomposition of a `Circulant` matrix \(A\).
   !    !!
   !    !! #### Syntax
   !    !!
   !    !! ```fortran
   !    !!    call svd(A, s, u, vt)
   !    !! ```
   !    !!
   !    !! #### Arguments
   !    !!
   !    !! `A`   :  Matrix of `Circulant` type.
   !    !!          It is an `intent(in)` argument.
   !    !!
   !    !! `s`   :  Rank-1 array `real` array returning the singular values of `A`.
   !    !!          It is an `intent(out)` argument.
   !    !!
   !    !! `u` (optional) :  Rank-2 array of the same kind as `A` returning the left singular
   !    !!                   vectors of `A` as columns. Its size should be `[n, n]`.
   !    !!                   It is an `intent(out)` argument.
   !    !!
   !    !! `vt (optional) :  Rank-2 array of the same kind as `A` returning the right singular
   !    !!                   vectors of `A` as rows. Its size should be `[n, n]`.
   !    !!                   It is an `intent(out)` argument.
   !    module subroutine svd_rdp(A, s, u, vt)
   !       !! Compute the singular value decomposition of a `Circulant` matrix.
   !       type(Circulant), intent(in) :: A
   !       !! Input matrix.
   !       real(dp), intent(out) :: s(:)
   !       !! Singular values in descending order.
   !       real(dp), optional, intent(out) :: u(:, :)
   !       !! Left singular vectors as columns.
   !       real(dp), optional, intent(out) :: vt(:, :)
   !       !! Right singular vectors as rows.
   !    end subroutine
   ! end interface

   !--------------------------------------------
   !-----     Eigenvalue Decomposition     -----
   !--------------------------------------------

   interface eigvals
      !! This interface overloads the `eigvals` interface from `stdlib_linalg`
      !! to compute the eigenvalues of a real-valued matrix \( A \) whose
      !! type is `Circulant`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  `real`-valued matrix of `Circulant` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda`  :  Vector of eigenvalues in increasing order.
      !!
      !! @note
      !! Eigenvalues of a circulant matrix can be efficiently computed using
      !! the Fast Fourier Transform of the generating vector `c`.
      !! @endnote
      module function eigvals_rdp(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real `Circulant`
         !! matrix.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         complex(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function
   end interface

   interface eig
      !! This interface overloads the `eig` interface from `stdlib_linalg` to
      !! compute the eigenvalues and eigenvectors of a real-valued matrix \(A\)
      !! whose type is `Circulant`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eig(A, lambda [, left] [, right])
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` : `real`-valued matrix of `Circulant`.
      !!          It is an `intent(in)` argument.
      !!
      !! - `lambda` :   Rank-1 `real` array returning the eigenvalues of `A`
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
      !! Eigenvalues of a circulant matrix can be efficiently computed using
      !! the Fast Fourier Transform of the generating vector `c`. Likewise,
      !! its eigenvectors are simply the corresponding Fourier modes.
      !! @endnote
      module subroutine eig_rdp(A, lambda, left, right)
         !! Utility function to compute the eigenvalues and eigenvectors of a
         !! `Circulant` matrix.
         type(Circulant), intent(in) :: A
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
      !! Convert a `Circulant` matrix to its dense representation.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = dense(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Circulant` type.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Rank-2 array representation of the matrix \( A \).
      module function dense_rdp(A) result(B)
         !! Utility function to convert a `Circulant` matrix to a rank-2 array.
         type(Circulant), intent(in) :: A
         !! Input diagonal matrix.
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function
   end interface

   interface transpose
      !! This interface overloads the Fortran `intrinsic` procedure to define
      !! the transpose operation of a `Circulant` matrix.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = transpose(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! - `A` :  Matrix of `Circulant`.
      !!          It is an `intent(in)` argument.
      !!
      !! - `B` :  Resulting transposed matrix. It is of the same type as `A`.
      pure module function transpose_rdp(A) result(B)
         !! Utility function to compute the transpose of a `Circulant` matrix.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         type(Circulant) :: B
         !! Transpose of the matrix.
      end function
   end interface

   interface size
      !! Utility function to return the size of a `Circulant` matrix along
      !! a given dimension.
      pure module function size_rdp(A, dim) result(arr_size)
         !! Utility function to return the size of `Circulant` matrix along a
         !! given dimension.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         integer(ilp), optional, intent(in) :: dim
         !! Queried dimension.
         integer(ilp) :: arr_size
         !! Size of the matrix along the dimension dim.
      end function
   end interface

   interface shape
      !! Utility function to return the shape of a `Circulant` matrix.
      pure module function shape_rdp(A) result(arr_shape)
         !! Utility function to get the shape of a `Circulant` matrix.
         type(Circulant), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: arr_shape(2)
         !! Shape of the matrix.
      end function
   end interface

   interface operator(*)
      pure module function scalar_multiplication_rdp(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a
         !! `Circulant` matrix.
         real(dp), intent(in) :: alpha
         type(Circulant), intent(in) :: A
         type(Circulant) :: B
      end function scalar_multiplication_rdp

      pure module function scalar_multiplication_bis_rdp(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a
         !! `Circulant` matrix.
         type(Circulant), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(Circulant) :: B
      end function scalar_multiplication_bis_rdp
   end interface
end module
