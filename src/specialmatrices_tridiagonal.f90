module SpecialMatrices_Tridiagonal
   use stdlib_kinds, only: wp => dp, int32
   implicit none
   private
   
   ! --> Linear Algebra.
   public :: transpose
   public :: matmul
   public :: solve
   ! public :: eig

   ! --> Utility functions.
   public :: dense

   !-----------------------------------------------------------
   !-----     Base types for bi/tri-diagonal matrices     -----
   !-----------------------------------------------------------

   type, public :: Diagonal
      !! Base type used to define a `Diagonal` matrix of size [n x n] with diagonal given by the
      !! rank-1 array `dv`.
      integer(int32) :: n
      !! Dimension of the matrix.
      real(wp), allocatable :: dv(:)
      !! Diagonal elements of the matrix.
   end type

   type, public :: Bidiagonal
      integer(int32) :: n
      !! Dimension of the matrix.
      real(wp), allocatable :: dv(:), ev(:)
      !! Diagonal elements
      character(len=1) :: which
      !! Whether it is lower- or upper-bidiagonal.
   end type

   type, public :: Tridiagonal
      integer(int32) :: n
      !! Dimension of the matrix.
      real(wp), allocatable :: d(:), du(:), dl(:)
      !! Tridiagonal elements.
   end type

   type, public :: SymTridiagonal
      integer(int32) :: n
      !! Dimension of the matrix.
      real(wp), allocatable :: dv(:), ev(:)
      !! Tridiagonal elements.
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
      pure module function initialize_diag(n) result(A)
         !! Utility function to construct a `Diagonal` matrix filled with zeros.
         integer(int32), intent(in) :: n
         !! Dimension of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function initialize_diag

      pure module function construct_diag(dv) result(A)
         !! Utility function to construct a `Diagonal` matrix from a rank-1 array.
         real(wp), intent(in) :: dv(:)
         !! Diagonal elements of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function construct_diag

      pure module function construct_constant_diag(d, n) result(A)
         !! Utility function to construct a `Diagonal` matrix with constant diagonal element.
         real(wp), intent(in) :: d
         !! Constant diagonal element of the matrix.
         integer(int32), intent(in) :: n
         !! Dimension of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function construct_constant_diag

      module function construct_dense_to_diag(A) result(B)
         !! Utility function to construct a `Diagonal` matrix from a rank-2 array.
         real(wp), intent(in) :: A(:, :)
         !! Dense [n x n] matrix from which to construct the `Diagonal` one.
         type(Diagonal) :: B
         !! Corresponding diagonal matrix.
      end function construct_dense_to_diag
   end interface





   interface Bidiagonal
      pure module function initialize_bidiag(n, which) result(A)
         integer(int32), intent(in) :: n
         character(len=1), optional, intent(in) :: which
         type(Bidiagonal) :: A
      end function initialize_bidiag

      pure module function construct_bidiag(dv, ev, which) result(A)
         real(wp), intent(in) :: dv(:), ev(:)
         character(len=1), optional, intent(in) :: which
         type(Bidiagonal) :: A
      end function construct_bidiag

      pure module function construct_constant_bidiag(d, e, n, which) result(A)
         real(wp), intent(in) :: d, e
         integer(int32), intent(in) :: n
         character(len=1), optional, intent(in) :: which
         type(Bidiagonal) :: A
      end function construct_constant_bidiag
   end interface




   interface Tridiagonal
      pure module function initialize_tridiag(n) result(A)
         integer(int32), intent(in) :: n
         type(Tridiagonal) :: A
      end function initialize_tridiag

      pure module function construct_tridiag(dl, d, du) result(A)
         real(wp), intent(in) :: dl(:), d(:), du(:)
         type(Tridiagonal) :: A
      end function construct_tridiag

      pure module function construct_constant_tridiag(l, d, u, n) result(A)
         real(wp), intent(in) :: l, d, u
         integer(int32), intent(in) :: n
         type(Tridiagonal) :: A
      end function construct_constant_tridiag
   end interface





   interface SymTridiagonal
      pure module function initialize_symtridiag(n) result(A)
         integer(int32), intent(in) :: n
         type(SymTridiagonal) :: A
      end function initialize_symtridiag

      pure module function construct_symtridiag(dv, ev) result(A)
         real(wp), intent(in) :: dv(:), ev(:)
         type(SymTridiagonal) :: A
      end function construct_symtridiag

      pure module function construct_constant_symtridiag(d, e, n) result(A)
         real(wp), intent(in) :: d, e
         integer(int32), intent(in) :: n
         type(SymTridiagonal) :: A
      end function construct_constant_symtridiag
   end interface

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   interface matmul
      !! This interface overloads the Fortran intrinsic `matmul` for the following types:
      !!
      !! - `Diagonal`
      !! - `Bidiagonal`
      !! - `Tridiagonal`
      !! - `SymTridiagonal`
      !!
      !! The intrinsic `matmul` is overloaded both for matrix-vector and matrix-matrix products.
      !! For a matrix-matrix product \( C = AB \), only the matrix \( A \) has to be of the ones
      !! of the types defined by `SpecialMatrices`. Both \( B \) and \( C \) need to be standard 
      !! Fortran rank-2 arrays. All the underlying functions are defined as `pure`.
      !!
      !! #### Syntax
      !!
      !! - For matrix-vector product with `A` being of a type defined by `SpecialMatrices` and
      !! `x` a standard rank-1 array:
      !! ```fortran
      !!    matmul(A, x)
      !! ```
      !!
      !! - For matrix-matrix product with `A` being of a type defined by `SpecialMatrices` and
      !! `B` a rank-2 array:
      !! ```fortran
      !!    matmul(A, B)
      !! ```
      pure module function diag_spmv(A, x) result(y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Diagonal` type and `x` and `y` are both rank-1 arrays.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(wp), intent(in) :: x(:)
         !! Input vector.
         real(wp) :: y(size(x))
         !! Output vector.
      end function diag_spmv

      pure module function diag_multi_spmv(A, X) result(Y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Diagonal` type and `X` and `Y` are both rank-2 arrays.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(wp), intent(in) :: X(:, :)
         !! Input matrix (rank-2 array).
         real(wp) :: Y(size(X, 1), size(X, 2))
         !! Output matrix (rank-2 array).
      end function diag_multi_spmv

      pure module function bidiag_spmv(A, x) result(y)
         type(Bidiagonal), intent(in) :: A
         real(wp), intent(in) :: x(:)
         real(wp) :: y(size(x))
      end function bidiag_spmv

      pure module function bidiag_multi_spmv(A, X) result(Y)
         type(Bidiagonal), intent(in) :: A
         real(wp), intent(in) :: X(:, :)
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function bidiag_multi_spmv

      pure module function tridiag_spmv(A, x) result(y)
         type(Tridiagonal), intent(in) :: A
         real(wp), intent(in) :: x(:)
         real(wp) :: y(size(x))
      end function tridiag_spmv

      pure module function tridiag_multi_spmv(A, X) result(Y)
         type(Tridiagonal), intent(in) :: A
         real(wp), intent(in) :: X(:, :)
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function tridiag_multi_spmv

      pure module function symtridiag_spmv(A, x) result(y)
         type(SymTridiagonal), intent(in) :: A
         real(wp), intent(in) :: x(:)
         real(wp) :: y(size(x))
      end function symtridiag_spmv

      pure module function symtridiag_multi_spmv(A, X) result(Y)
         type(SymTridiagonal), intent(in) :: A
         real(wp), intent(in) :: X(:, :)
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function symtridiag_multi_spmv
   end interface

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   interface solve
      !! This interface overloads the `solve` interface from `stdlib_linalg` for solving a linear
      !! system \( Ax = b \) where \( A \) is of one of the types provided by `SpecialMatrices`.
      !! It also enables to solve a linear system with multiple right-hand sides.
      !!
      !! #### Syntax
      !!
      !! To solve a system with \( A \) being of one of the types defined by `SpecialMatrices`:
      !! 
      !! ```fortran
      !!    y = solve(A, x)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `x`   :  Rank-1 or rank-2 array defining the right-hand side(s). It is an `intent(in)`
      !!          argument.
      !!
      !! `y`   :  Solution of the linear system. It has the same type and shape as `x`.
      pure module function diag_solve(A, b) result(x)
         !! Utility function to solve the linear system \( A x = b \) where \( A \) is of
         !! `Diagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
         !! with the same type and dimension as `b`.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(wp), intent(in) :: b(:)
         !! Right-hande side vector.
         real(wp) :: x(size(b))
         !! Solution vector.
      end function diag_solve

      pure module function diag_multi_solve(A, B) result(X)
         !! Utility function to solve a linear system with multiple right-hand sides where
         !! \( A \) is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a
         !! rank-2 array with the same type and dimensions as `B`.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(wp), intent(in) :: B(:, :)
         !! Right-hande side vectors.
         real(wp) :: X(size(B, 1), size(B, 2))
         !! Solution vectors.
      end function diag_multi_solve


      ! Bidiagonal matrix solve.
      pure module function bidiag_solve(A, b) result(x)
         type(Bidiagonal), intent(in) :: A
         real(wp), intent(in) :: b(:)
         real(wp) :: x(size(b))
      end function bidiag_solve

      pure module function bidiag_multi_solve(A, B) result(X)
         type(Bidiagonal), intent(in) :: A
         real(wp), intent(in) :: B(:, :)
         real(wp) :: X(size(B, 1), size(B, 2))
      end function bidiag_multi_solve

      
      ! Tridiagonal matrix solve.
      pure module function tridiag_solve(A, b) result(x)
         type(Tridiagonal), intent(in) :: A
         real(wp), intent(in) :: b(:)
         real(wp) :: x(size(b))
      end function tridiag_solve

      pure module function tridiag_multi_solve(A, B) result(X)
         type(Tridiagonal), intent(in) :: A
         real(wp), intent(in) :: B(:, :)
         real(wp) :: X(size(B, 1), size(B, 2))
      end function tridiag_multi_solve

      
      ! Symmetric Tridiagonal matrix solve.
      pure module function symtridiag_solve(A, b) result(x)
         type(SymTridiagonal), intent(in) :: A
         real(wp), intent(in) :: b(:)
         real(wp) :: x(size(b))
      end function symtridiag_solve

      pure module function symtridiag_multi_solve(A, B) result(X)
         type(SymTridiagonal), intent(in) :: A
         real(wp), intent(in) :: B(:, :)
         real(wp) :: X(size(B, 1), size(B, 2))
      end function symtridiag_multi_solve
   end interface




   interface eig
      !! This interface overloads the `eig` interface from `stdlib_linalg` for solving a
      !! generalized eigenvalue problem \( A x = \lambda x \) where \( A \) is of one of the types
      !! provided by `SpecialMatrices`. Whenever possible, eigensolvers specialized for the
      !! particular matrix structure are being used. If not, \( A \) is converted to a standard
      !! rank-2 array and the eigenpairs are computed using the `eig` from `stdlib_linalg`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eig(A, lambda, vectors)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type. It
      !!          is an `intent(in)` argument.
      !!
      !! `lambda` :  `complex` rank-1 array containing the eigenvalues of \( A \). It is an
      !!             `intent(out)` argument.
      !!
      !! `vectors`   :  `complex` rank-2 array containing th eigenvectors of \( A \). It is an
      !!                `intent(out)` argument.
      pure module subroutine diag_eig(A, lambda, vectors)
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(wp), intent(out) :: lambda(A%n)
         !! Eigenvalues.
         real(wp), intent(out) :: vectors(A%n, A%n)
         !! Eigenvectors.
      end subroutine diag_eig
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
      pure module function diag_to_dense(A) result(B)
         !! Utility function to convert a `Diagonal` matrix to a regular rank-2 array.
         type(Diagonal), intent(in) :: A
         !! Input diagonal matrix.
         real(wp) :: B(A%n, A%n)
         !! Output dense rank-2 array.
      end function diag_to_dense

      pure module function bidiag_to_dense(A) result(B)
         type(Bidiagonal), intent(in) :: A
         real(wp) :: B(A%n, A%n)
      end function bidiag_to_dense

      pure module function tridiag_to_dense(A) result(B)
         type(Tridiagonal), intent(in) :: A
         real(wp) :: B(A%n, A%n)
      end function tridiag_to_dense

      pure module function symtridiag_to_dense(A) result(B)
         type(SymTridiagonal), intent(in) :: A
         real(wp) :: B(A%n, A%n)
      end function symtridiag_to_dense
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
      pure module function diag_transpose(A) result(B)
         !! Utility function to compute the transpose of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input Diagonal matrix.
         type(Diagonal) :: B
         !! Tranpose of the original diagonal matrix.
      end function diag_transpose

      pure module function bidiag_transpose(A) result(B)
         type(Bidiagonal), intent(in) :: A
         type(Bidiagonal) :: B
      end function bidiag_transpose

      pure module function tridiag_transpose(A) result(B)
         type(Tridiagonal), intent(in) :: A
         type(Tridiagonal) :: B
      end function tridiag_transpose

      pure module function symtridiag_transpose(A) result(B)
         type(SymTridiagonal), intent(in) :: A
         type(SymTridiagonal) :: B
      end function symtridiag_transpose
   end interface

contains

end module SpecialMatrices_Tridiagonal
