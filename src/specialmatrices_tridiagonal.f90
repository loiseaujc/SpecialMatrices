module SpecialMatrices_Tridiagonal
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling
   implicit none
   private

   ! --> Linear Algebra.
   public :: transpose
   public :: det
   public :: trace
   public :: inv
   public :: matmul, spmv_ip
   public :: solve, solve_ip
   public :: svd, svdvals
   public :: eigh, eigvalsh

   ! --> Utility functions.
   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   !-----------------------------------------------------------
   !-----     Base types for bi/tri-diagonal matrices     -----
   !-----------------------------------------------------------

   type, public :: Diagonal
      !! Base type used to define a `Diagonal` matrix of size [n x n] with diagonal given by the
      !! rank-1 array `dv`.
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: dv(:)
      !! Diagonal elements of the matrix.
   end type

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

   type, public :: SymTridiagonal
      !! Base type used to define a `SymTridiagonal` matrix of size [n x n] with diagonal elements
      !! given by the rank-1 arrays `dv` (diagonal) and `ev` (sub- and super-diagonal).
      private
      integer(ilp) :: n
      !! Dimension of the matrix.
      real(dp), allocatable :: dv(:), ev(:)
      !! Tridiagonal elements.
      logical :: isposdef
      !! Whether `A` is symmetric positive definite or not.
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
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function initialize_diag

      pure module function construct_diag(dv) result(A)
         !! Utility function to construct a `Diagonal` matrix from a rank-1 array.
         real(dp), intent(in) :: dv(:)
         !! Diagonal elements of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function construct_diag

      pure module function construct_constant_diag(d, n) result(A)
         !! Utility function to construct a `Diagonal` matrix with constant diagonal element.
         real(dp), intent(in) :: d
         !! Constant diagonal element of the matrix.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(Diagonal) :: A
         !! Corresponding diagonal matrix.
      end function construct_constant_diag

      module function construct_dense_to_diag(A) result(B)
         !! Utility function to construct a `Diagonal` matrix from a rank-2 array.
         real(dp), intent(in) :: A(:, :)
         !! Dense [n x n] matrix from which to construct the `Diagonal` one.
         type(Diagonal) :: B
         !! Corresponding diagonal matrix.
      end function construct_dense_to_diag
   end interface

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

   interface SymTridiagonal
      !! This interface provides different methods to construct a `SymTridiagonal` matrix.
      !! Only `double precision` is supported currently. Only the non-zero elements of
      !! \( A \) are stored, i.e.
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
      !! - If \( A \) is known to be symmetric positive definite, it can be constructed as
      !!
      !! ```fortran
      !!    A = SymTridiagonal(dv, ev, ifposdef=.true.)
      !! ```
      pure module function initialize_symtridiag(n) result(A)
         !! Utility function to create a `SymTridiagonal` matrix filled with zeros.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         type(SymTridiagonal) :: A
         !! Corresponding symmetric tridiagonal matrix.
      end function initialize_symtridiag

      pure module function construct_symtridiag(dv, ev, isposdef) result(A)
         !! Utility function to create a `SymTridiagonal` matrix from rank-1 arrays.
         real(dp), intent(in) :: dv(:), ev(:)
         !! Diagonal elements of the matrix.
         logical, optional, intent(in) :: isposdef
         !! Whether `A` is symmetric positive definite or not.
         type(SymTridiagonal) :: A
         !! Corresponding symmetric tridiagonal matrix.
      end function construct_symtridiag

      pure module function construct_constant_symtridiag(d, e, n, isposdef) result(A)
         !! Utility function to create a `SymTridiagonal` matrix with constant diagonal elements.
         real(dp), intent(in) :: d, e
         !! Constant diagonal elements of the matrix.
         integer(ilp), intent(in) :: n
         !! Dimension of the matrix.
         logical, optional, intent(in) :: isposdef
         !! Whether `A` is symmetric positive definite or not.
         type(SymTridiagonal) :: A
         !! Corresponding symmetric tridiagonal matrix.
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
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function diag_spmv

      pure module function diag_multi_spmv(A, X) result(Y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Diagonal` type and `X` and `Y` are both rank-2 arrays.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: X(:, :)
         !! Input matrix (rank-2 array).
         real(dp), allocatable :: Y(:, :)
         !! Output matrix (rank-2 array).
      end function diag_multi_spmv

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

      pure module function symtridiag_spmv(A, x) result(y)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `SymTridiagonal` type and `x` and `y` are both rank-1 arrays.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
         real(dp), allocatable :: y(:)
         !! Output vector.
      end function symtridiag_spmv

      pure module function symtridiag_multi_spmv(A, X) result(Y)
         !! Utility function to compute the matrix-matrix product \( Y = AX \) where \( A \)
         !! is of `SymTridiagonal` type and `X` and `Y` are both rank-2 arrays.
         type(SymTridiagonal), intent(in) :: A
         real(dp), intent(in) :: X(:, :)
         real(dp), allocatable :: Y(:, :)
      end function symtridiag_multi_spmv
   end interface

   interface spmv_ip
      !! This interface provides methods for in-place matrix-vector products for the following
      !! types:
      !!
      !! - `Diagonal`
      !! - `Bidiagonal`
      !! - `Tridiagonal`
      !! - `SymTridiagonal`
      !!
      !! For a matrix-matrix product \( C = AB \), only the matrix \( A \) has to be of the ones
      !! of the types defined by `SpecialMatrices`. Both \( B \) and \( C \) need to be standard
      !! Fortran rank-2 arrays. All the underlying functions perform the computations in-place,
      !! the array `y`/`C` will be overwritten with the result.
      !!
      !! #### Syntax
      !!
      !! - For matrix-vector product with `A` being of a type defined by `SpecialMatrices` and
      !! `x` a standard rank-1 array:
      !! ```fortran
      !!    call spmv_ip(y, A, x)
      !! ```
      !!
      !! - For matrix-matrix product with `A` being of a type defined by `SpecialMatrices` and
      !! `B` a rank-2 array:
      !! ```fortran
      !!    call spmv_ip(C, A, B)
      !! ```
      module subroutine diag_spmv_ip(y, A, x)
         !! Utility function to compute the matrix-vector product \( y = Ax \) where \( A \)
         !! is of `Diagonal` type and `x` and `y` are both rank-1 arrays. Note that this
         !! function performs this product in-place, i.e. `y` needs to be pre-allocated and
         !! will be modified by the call.
         real(dp), intent(out) :: y(:)
         !! Output vector.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: x(:)
         !! Input vector.
      end subroutine diag_spmv_ip

      module subroutine diag_multi_spmv_ip(Y, A, X)
         !! Utility function to compute the matrix-matrix product \( Y = AX \) where \( A \)
         !! is of `Diagonal` type and `X` and `Y` are both rank-2 arrays.
         real(dp), intent(out) :: Y(:, :)
         !! Output vectors.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), intent(in) :: X(:, :)
         !! Input vectors.
      end subroutine diag_multi_spmv_ip
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
      pure module function diag_solve(A, b) result(x)
         !! Utility function to solve the linear system \( A x = b \) where \( A \) is of
         !! `Diagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1 array
         !! with the same type and dimension as `b`.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hande side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function diag_solve

      pure module function diag_multi_solve(A, B) result(X)
         !! Utility function to solve a linear system with multiple right-hand sides where
         !! \( A \) is of `Diagonal` type and `B` a rank-2 array. The solution `X` is also a
         !! rank-2 array with the same type and dimensions as `B`.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hande side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function diag_multi_solve

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

      ! Symmetric Tridiagonal matrix solve.
      pure module function symtridiag_solve(A, b) result(x)
         !! Utility function to solve the linear system \( Ax = b \) where \( A \) is of
         !! `SymTridiagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1
         !! array with the same type and dimension as `b`.
         type(SymTridiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hande side vector.
         real(dp), allocatable :: x(:)
         !! Solution vector.
      end function symtridiag_solve

      pure module function symtridiag_multi_solve(A, B) result(X)
         !! Utility function to solve a linear system with multiple right-hand side vectors
         !! where \( A \) is of `SymTridiagonal` type and `B` a rank-2 array. The solution `X`
         !! is also a rank-2 array with the same type and dimensions as `B`.
         type(SymTridiagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: B(:, :)
         !! Right-hand side vectors.
         real(dp), allocatable :: X(:, :)
         !! Solution vectors.
      end function symtridiag_multi_solve
   end interface

   interface solve_ip
      !! This interface provides methods for solving *in-place* a linear
      !! system \( Ax = b \) where \( A \) is of one of the types provided by `SpecialMatrices`.
      !! It also enables to solve a linear system with multiple right-hand sides.
      !!
      !! #### Syntax
      !!
      !! To solve a system with \( A \) being of one of the types defined by `SpecialMatrices`:
      !!
      !! ```fortran
      !!    call solve_ip(x, A, b)
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
      module subroutine diag_solve_ip(x, A, b)
         !! Utility function to solve the linear system \( Ax = b \) where \( A \) is of
         !! `Diagonal` type and `b` a rank-1 array. The solution `x` is also a rank-1
         !! array with the same type and dimension as `b`. Computation is done in-place, i.e.
         !! the array `x` will be overwritten with the solution.
         real(dp), intent(out) :: x(:)
         !! Solution vector.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:)
         !! Right-hand side vector.
      end subroutine diag_solve_ip

      module subroutine diag_multi_solve_ip(x, A, b)
         !! Utility function to solve the linear system \( Ax = b \) where \( A \) is of
         !! `Diagonal` type and `B` a rank-2 array. The solution `x` is also a rank-2
         !! array with the same type and dimension as `B`. Computation is done in-place, i.e.
         !! the array `x` will be overwritten with the solution.
         real(dp), intent(out) :: x(:, :)
         !! Solution vector.
         type(Diagonal), intent(in) :: A
         !! Coefficient matrix.
         real(dp), intent(in) :: b(:, :)
         !! Right-hand side vector.
      end subroutine diag_multi_solve_ip
   end interface

   interface det
      !! This interface overloads the `det` interface from `stdlib_linag` to compute the
      !! determinant \(\det(A)\) where \(A\) is of one of the types provided by `SpecialMatrices`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    d = det(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is in an `intent(in)` argument.
      !!
      !! `determinant`  :  Determinant of the matrix.
      pure module function diag_det(A) result(determinant)
         !! Utility function to compute the determinant of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: determinant
         !! Determinant of the matrix.
      end function diag_det
   end interface

   interface trace
      !! This interface overloads the `trace` interface from `stdlib_linalg` to compute the trace
      !! of a matrix \( A \) whose type is one of the types provided by `SpecialMatrices`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    tr = trace(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `tr`  :  Trace of the matrix.
      pure module function diag_trace(A) result(tr)
         !! Utility function to compute the trace of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp) :: tr
         !! Trace of the matrix.
      end function diag_trace
   end interface

   interface inv
      !! This interface overloads the `inv` interface from `stdlib_linalg` to compute the inverse
      !! of a matrix \( A \) whose type is one of the types provided by `SpecialMatrices`. Note
      !! that the inverse is returned as a standard Fortran rank-2 array.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    B = inv(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `B`   :  Inverse of the matrix \(A\). It is a standard Fortran rank-2 array.
      pure module function diag_inv(A) result(B)
         !! Utility function to compute the inverse of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: B(:, :)
         !! Inverse of the matrix.
      end function diag_inv
   end interface

   interface svdvals
      !! This interface overloads the `svdvals` interface from `stdlib_linalg` to compute the
      !! singular values of a matrix \( A \) whose type is one of the types provided by
      !! `SpecialMatrices`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    s = svdvals(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `s`   :  Vector of singular values sorted in decreasing order.
      pure module function diag_svdvals(A) result(s)
         !! Utility function to compute the singular values of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: s(:)
         !! Singular values.
      end function diag_svdvals
   end interface

   interface svd
      !! This interface overloads the `svd` interface from `stdlib_linalg` to compute the
      !! the singular value decomposition of a matrix \( A \) whose type is provided by
      !! `SpecialMatrices`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call svd(A, s, u, vt)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  Matrix of `Diagonal`, `Bidiagonal`, `Tridiagonal` or `SymTridiagonal` type.
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
      module subroutine diag_svd(A, u, s, vt)
         !! Utility function to compute the singular value decomposition of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out) :: s(:)
         !! Singular values.
         real(dp), allocatable, optional, intent(out) :: u(:, :)
         !! Left singular vectors as columns.
         real(dp), allocatable, optional, intent(out) :: vt(:, :)
         !! Right singular vectors as rows.
      end subroutine diag_svd
   end interface

   interface eigvalsh
      !! This interface overloads the `eigvalsh` interface from `stdlib_linalg` to compute the
      !! eigenvalues of a real-valued matrix \( A \) whose type is `Diagonal` or `SymTridiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    lambda = eigvalsh(A)
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   :  `real`-valued matrix of `Diagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `lambda` :  Vector of eigenvalues in increasing order.
      module function diag_eigvalsh(A) result(lambda)
         !! Utility function to compute the eigenvalues of a real `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable :: lambda(:)
         !! Eigenvalues.
      end function diag_eigvalsh
   end interface

   interface eigh
      !! This interface overloads the `eigh` interface from `stdlib_linalg` to compute the
      !! eigenvalues and eigenvectors of a real-valued matrix \(A\) whose type is `Diagonal`
      !! or `SymTridiagonal`.
      !!
      !! #### Syntax
      !!
      !! ```fortran
      !!    call eigh(A, lambda [, vectors])
      !! ```
      !!
      !! #### Arguments
      !!
      !! `A`   : `real`-valued matrix of `Diagonal` or `SymTridiagonal` type.
      !!          It is an `intent(in)` argument.
      !!
      !! `lambda` :  Rank-1 `real` array returning the eigenvalues of `A` in increasing order.
      !!             It is an `intent(out)` argument.
      !!
      !! `vectors` (optional) :  Rank-2 array of the same kind as `A` returning the eigenvectors
      !!                         of `A`. It is an `intent(out)` argument.
      module subroutine diag_eigh(A, lambda, vectors)
         !! Utility function to compute the eigenvalues and eigenvectors of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         real(dp), allocatable, intent(out) :: lambda(:)
         !! Eigenvalues.
         real(dp), allocatable, optional, intent(out) :: vectors(:, :)
         !! Eigenvectors.
      end subroutine diag_eigh
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
         real(dp), allocatable :: B(:, :)
         !! Output dense rank-2 array.
      end function diag_to_dense

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

      pure module function symtridiag_to_dense(A) result(B)
         !! Utility function to convert a `SymTridiagonal` matrix to a regular rank-2 array.
         type(SymTridiagonal), intent(in) :: A
         !! Input symmetric tridiagonal matrix.
         real(dp) :: B(A%n, A%n)
         !! Output dense rank-2 array.
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
         !! The output matrix is also of `Diagonal` type.
         type(Diagonal), intent(in) :: A
         !! Input Diagonal matrix.
         type(Diagonal) :: B
         !! Tranpose of the original diagonal matrix.
      end function diag_transpose

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

      pure module function symtridiag_transpose(A) result(B)
         !! Utility function to compute the transpose of a `SymTridiagonal` matrix.
         !! The output matrix is also of `SymTridiagonal` type.
         type(SymTridiagonal), intent(in) :: A
         !! Input symmetric tridiagonal matrix.
         type(SymTridiagonal) :: B
         !! Transpose of the original matrix.
      end function symtridiag_transpose
   end interface

   interface size
      pure module function diag_size(A, dim) result(arr_size)
         !! Utility function to return the size of a `Diagonal` matrix along a given dimension.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp), intent(in) :: dim
         !! Dimension whose size needs to be known.
         integer(ilp) :: arr_size
      end function diag_size

      pure module function symtridiag_size(A, dim) result(arr_size)
         !! Utility function to return the size of a `Diagonal` matrix along a given dimension.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp), intent(in) :: dim
         !! Dimension whose size needs to be known.
         integer(ilp) :: arr_size
      end function symtridiag_size
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
      pure module function diag_shape(A) result(shape)
         !! Utility function to get the shape of a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: shape(2)
         !! Shape of the matrix.
      end function diag_shape

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

      pure module function symtridiag_shape(A) result(shape)
         !! Utility function to get the shape of a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         !! Input matrix.
         integer(ilp) :: shape(2)
         !! Shape of the matrix.
      end function symtridiag_shape
   end interface

   interface operator(*)
      pure module function diag_scalar_multiplication(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a `Diagonal` matrix.
         real(dp), intent(in) :: alpha
         type(Diagonal), intent(in) :: A
         type(Diagonal) :: B
      end function diag_scalar_multiplication

      pure module function diag_scalar_multiplication_bis(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a `Diagonal` matrix.
         type(Diagonal), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(Diagonal) :: B
      end function diag_scalar_multiplication_bis

      pure module function symtridiag_scalar_multiplication(alpha, A) result(B)
         !! Utility function to perform a scalar multiplication with a `SymTridiagonal` matrix.
         real(dp), intent(in) :: alpha
         type(SymTridiagonal), intent(in) :: A
         type(SymTridiagonal) :: B
      end function symtridiag_scalar_multiplication

      pure module function symtridiag_scalar_multiplication_bis(A, alpha) result(B)
         !! Utility function to perform a scalar multiplication with a `SymTridiagonal` matrix.
         type(SymTridiagonal), intent(in) :: A
         real(dp), intent(in) :: alpha
         type(SymTridiagonal) :: B
      end function symtridiag_scalar_multiplication_bis
   end interface

contains

end module SpecialMatrices_Tridiagonal
