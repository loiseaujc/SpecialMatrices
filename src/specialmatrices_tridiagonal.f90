module SpecialMatrices_Tridiagonal
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling
   implicit none
   private

   ! --> Linear Algebra.
   public :: transpose
   public :: matmul
   public :: solve

   ! --> Utility functions.
   public :: dense

   !-----------------------------------------------------------
   !-----     Base types for bi/tri-diagonal matrices     -----
   !-----------------------------------------------------------

   type, public :: Bidiagonal
      !! Base type used to define a `Bidiagonal` matrix of size [n x n] with diagonals given by
      !! the rank-1 arrays `dv` and `ev`. The character `which` determines whether `ev` defines the
      !! sub-diagonal (`which = "L"`, default) or the super-diagonal (`which = "U"`).
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: dv(:), ev(:)
      !! Diagonal elements
      character(len=1) :: which
      !! Whether it is lower- or upper-bidiagonal.
   end type

   type, public :: Tridiagonal
      !! Base type used to define a `Tridiagonal` matrix of size [n x n] with diagonal elements
      !! given by the rank-1 arrays `dl` (sub-diagonal), `d` (diagonal) and `du` (super-diagonal).
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: d(:), du(:), dl(:)
      !! Tridiagonal elements.
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Bidiagonal
      !! This interface provides different methods to construct a `Bidiagonal` matrix.
      !! Only `double precision` is supported currently. Only the non-zero elements of
      !! \( A \) are stored, i.e.
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       d_1   \\
      !!       e_1  &  d_2   \\
      !!             &  \ddots   &  \ddots   \\
      !!             &           &  e_{n-2} &  d_{n-1}  \\
      !!             &           &           &  e_{n-1} &  d_n
      !!    \end{bmatrix}.
      !! \]
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
      !!    real(dp), allocatable :: dv(:), ev(:)
      !!    type(Bidiagonal) :: A
      !!    integer :: i
      !!
      !!    ev = [(i, i=1, n)]; dv = [(2*i, i=1, n)]
      !!    A = Bidiagonal(dv, ev, which="L")
      !! ```
      !!
      !! - Construct a `Bidiagonal` matrix with constant diagonals:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), parameter :: e = 1.0_dp, d = 2.0_dp
      !!    type(Bidiagonal) :: A
      !!
      !!    A = Bidiagonal(d, e, n, which="L")
      !! ```
      pure module function initialize_bidiag(n, which) result(A)
         !! Utility function to construct a `Bidiagonal` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         character(len=1), optional, intent(in) :: which
         !! Whether `A` has a sub- or super-diagonal.
         type(Bidiagonal) :: A
         !! Corresponding bidiagonal matrix.
      end function initialize_bidiag

      pure module function construct_bidiag(dv, ev, which) result(A)
         !! Utility function to construct a `Bidiagonal` matrix from rank-1 arrays.
         real(dp), intent(in) :: dv(:), ev(:)
         !! Diagonal elements of the matrix.
         character(len=1), optional, intent(in) :: which
         !! Whether `A` has a sub- or super-diagonal.
         type(Bidiagonal) :: A
         !! Corresponding bidiagonal matrix.
      end function construct_bidiag

      pure module function construct_constant_bidiag(d, e, n, which) result(A)
         !! Utility function to construct a `Bidiagonal` matrix with constant elements.
         real(dp), intent(in) :: d, e
         !! Constant diagonal elements.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         character(len=1), optional, intent(in) :: which
         !! Whether `A` has a sub- or super-diagonal.
         type(Bidiagonal) :: A
         !! Corresponding bidiagonal matrix.
      end function construct_constant_bidiag
   end interface

   interface Tridiagonal
      !! This interface provides different methods to construct a `Tridiagonal` matrix.
      !! Only `double precision` is supported currently. Only the non-zero elements of
      !! \( A \) are stored, i.e.
      !!
      !! \[
      !!    A
      !!    =
      !!    \begin{bmatrix}
      !!       d_1   &  du_1  \\
      !!       dl_1  &  d_2      &  du_2  \\
      !!             &  \ddots   &  \ddots   &  \ddots   \\
      !!             &           &  dl_{n-2} &  d_{n-1}  &  du_{n-1} \\
      !!             &           &           &  dl_{n-1} &  d_n
      !!    \end{bmatrix}.
      !! \]
      !!
      !! #### Syntax
      !!
      !! - Construct a `Tridiagonal` matrix filled with zeros:
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    type(Tridiagonal) :: A
      !!
      !!    A = Tridiagonal(n)
      !! ```
      !!
      !! - Construct a `Tridiagonal` matrix from rank-1 arrays:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), allocatable :: dl(:), d(:), du(:)
      !!    type(Tridiagonal) :: A
      !!    integer :: i
      !!
      !!    dl = [(i, i=1, n)]; d = [(2*i, i=1, n)]; dl = [(3*i, i=1, n)]
      !!    A = Tridiagonal(dl, d, du)
      !! ```
      !!
      !! - Construct a `Tridiagonal` matrix with constant diagonals:
      !!
      !! ```fortran
      !!    integer, parameter :: n
      !!    real(dp), parameter :: dl = 1.0_dp, d = 2.0_dp, du = 3.0_dp
      !!    type(Tridiagonal) :: A
      !!
      !!    A = Tridiagonal(dl, d, du, n)
      !! ```
      !!
      !! - Construct a `Tridiagonal` matrix from a standard Fortran rank-2 array:
      !!
      !! ```fortran
      !!    integer, parameter :: n = 100
      !!    real(dp) :: B(n, n)
      !!    type(Tridiagonal) :: A
      !!
      !!    call random_number(B); A = Tridiagonal(B)
      !! ```
      pure module function initialize_tridiag(n) result(A)
         !! Utility function to construct a `Tridiagonal` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Tridiagonal) :: A
         !! Corresponding tridiagonal matrix.
      end function initialize_tridiag

      pure module function construct_tridiag(dl, d, du) result(A)
         !! Utility function to construct a `Tridiagonal` matrix from a set of rank-1 arrays.
         real(dp), intent(in) :: dl(:), d(:), du(:)
         !! Diagonal elements of the matrix.
         type(Tridiagonal) :: A
         !! Corresponding tridiagonal matrix.
      end function construct_tridiag

      pure module function construct_constant_tridiag(l, d, u, n) result(A)
         !! Utility function to construct a `Tridiagonal` matrix with constant elements.
         real(dp), intent(in) :: l, d, u
         !! Constant diagonal elements.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Tridiagonal) :: A
         !! Corresponding tridiagonal matrix.
      end function construct_constant_tridiag

      module function construct_dense_to_tridiag(A) result(B)
         !! Utility function to construct a `Tridiagonal` matrix from a rank-2 array.
         real(dp), intent(in) :: A(:, :)
         !! Dense [n x n] matrix from which to construct the `Tridiagonal` one.
         type(Tridiagonal) :: B
         !! Corresponding tridiagonal matrix.
      end function construct_dense_to_tridiag
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
      pure module function bidiag_spmv(A, x) result(y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Bidiagonal` type and `x` and `y` are both rank-1 arrays.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function bidiag_spmv

      pure module function bidiag_multi_spmv(A, X) result(Y)
         !! Utility function to compute the matrix-matrix product \( Y = AX \) where \( A \)
         !! is of `Bidiagonal` type and `X` and `Y` are both rank-2 arrays.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: X(:, :)
         !! Input matrix (rank-2 array).
         real(dp), allocatable :: Y(:, :)
         !! Output matrix (rank-2 array).
      end function bidiag_multi_spmv

      pure module function tridiag_spmv(A, x) result(y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Tridiagonal` type and `x` and `y` are both rank-1 arrays.
         type(Tridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function tridiag_spmv

      pure module function tridiag_multi_spmv(A, X) result(Y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Tridiagonal` type and `X` and `Y` are both rank-2 arrays.
         type(Tridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: X(:, :)
         !! Input matrix (rank-2 array).
         real(dp), allocatable :: Y(:, :)
         !! Output matrix (rank-2 array).
      end function tridiag_multi_spmv
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
      !!    x = solve(A, b)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `b`   :  Rank-1 or rank-2 array defining the right-hand side(s). It is an `intent(in)`
      !!          argument.
      !!
      !! `x`   :  Solution of the linear system. It has the same type and shape as `b`.
      ! Bidiagonal matrix solve.
      pure module function bidiag_solve(A, b) result(x)
         !! Utility function to solve a linear system \( Ax = b \) where \( A \) is of
         !! `Bidiagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
         !! with the same type and dimension as `b`.
         type(Bidiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function bidiag_solve

      pure module function bidiag_multi_solve(A, B) result(X)
         !! Utility function to solve a linear system with multiple right-hand sides where
         !! \( A \) is of `Bidiagonal` type and `B` a rank-2 array. The solution `X` is also
         !! a rank-2 array with the same type and dimensions as `B`.
         type(Bidiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function bidiag_multi_solve

      ! Tridiagonal matrix solve.
      pure module function tridiag_solve(A, b) result(x)
         !! Utility function to solve the linear system \( Ax = b \) where \( A \) is of
         !! `Tridiagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1
         !! array with the same type and dimension as `b`.
         type(Tridiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function tridiag_solve

      pure module function tridiag_multi_solve(A, B) result(X)
         !! Utility function to solve a linear system with multiple right-hand sides where
         !! \( A \) is of `Tridiagonal` type and `B` a rank-2 array. The solution `X` is also
         !! a rank-2 array with the same type and dimensions as `B`.
         type(Tridiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function tridiag_multi_solve
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
      pure module function bidiag_to_dense(A) result(B)
         !! Utility function to convert a `Bidiagonal` matrix to a regular rank-2 array.
         type(Bidiagonal), intent(in) :: A
         !! Input bidiagonal matrix.
         real(dp) :: B(A%n, A%n)
         !! Output dense rank-2 array.
      end function bidiag_to_dense

      pure module function tridiag_to_dense(A) result(B)
         !! Utility function to convert a `Tridiagonal` matrix to a regular rank-2 array.
         type(Tridiagonal), intent(in) :: A
         !! Input tridiagonal matrix.
         real(dp) :: B(A%n, A%n)
         !! Output dense rank-2 array.
      end function tridiag_to_dense
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
      pure module function bidiag_transpose(A) result(B)
         !! Utility function to compute the transpose of a `Bidiagonal` matrix.
         !! The output matrix is also of `Bidiagonal` type.
         type(Bidiagonal), intent(in) :: A
         !! Input bidiagonal matrix.
         type(Bidiagonal) :: B
         !! Transpose of the original diagonal matrix.
      end function bidiag_transpose

      pure module function tridiag_transpose(A) result(B)
         !! Utility function to compute the tranpose of a `Tridiagonal` matrix.
         !! The output matrix is also of `Tridiagonal` type.
         type(Tridiagonal), intent(in) :: A
         !! Input tridiagonal matrix.
         type(Tridiagonal) :: B
         !! Transpose of the original tridiagonal matrix.
      end function tridiag_transpose
   end interface

   interface shape
      !! This interface provides methods to access the shape of a matrix \( A \) of one of
      !! the types defined by `SpecialMatrices`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    shape(A)
      !! ```
      pure module function bidiag_shape(A) result(shape)
         !! Utility function to get the shape of a `Bidiagonal` matrix.
         type(Bidiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: shape(2)
         !! Shape of the matrix.
      end function bidiag_shape

      pure module function tridiag_shape(A) result(shape)
         !! Utility function to get the shape of a `Tridiagonal` matrix.
         type(Tridiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: shape(2)
         !! Shape of the matrix.
      end function tridiag_shape
   end interface

contains

end module SpecialMatrices_Tridiagonal
