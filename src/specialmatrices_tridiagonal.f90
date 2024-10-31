module SpecialMatrices_Tridiagonal
   use stdlib_kinds, only: wp => dp, int32
   implicit none
   private
   
   ! --> Linear Algebra.
   public :: matmul
   public :: solve
   public :: transpose

   ! --> Utility functions.
   public :: dense

   !-----------------------------------------------------------
   !-----     Base types for bi/tri-diagonal matrices     -----
   !-----------------------------------------------------------

   type, public :: Diagonal
      ! Dimension of the matrix.
      integer(int32) :: n
      ! Diagonal elements.
      real(wp), allocatable :: dv(:)
   end type

   type, public :: Bidiagonal
      ! Dimension of the matrix.
      integer(int32) :: n
      ! Diagonal elements.
      real(wp), allocatable :: dv(:), ev(:)
      ! Whether it is lower- or upper-bidiagonal.
      character(len=1) :: which
   end type

   type, public :: Tridiagonal
      ! Dimension of the matrix.
      integer(int32) :: n
      ! Tridiagonal elements.
      real(wp), allocatable :: d(:), du(:), dl(:)
   end type

   type, public :: SymTridiagonal
      ! Dimension of the matrix.
      integer(int32) :: n
      ! Tridiagonal elements.
      real(wp), allocatable :: dv(:), ev(:)
   end type

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   interface Diagonal
      pure module function initialize_diag(n) result(A)
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Output matrix.
         type(Diagonal) :: A
      end function initialize_diag

      pure module function construct_diag(dv) result(A)
         ! Diagonal elements.
         real(wp), intent(in) :: dv(:)
         ! Output matrix.
         type(Diagonal) :: A
      end function construct_diag

      pure module function construct_constant_diag(d, n) result(A)
         ! Diagonal elements.
         real(wp), intent(in) :: d
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Output matrix.
         type(Diagonal) :: A
      end function construct_constant_diag

      module function construct_dense_to_diag(A) result(B)
         ! Input dense matrix.
         real(wp), intent(in) :: A(:, :)
         ! Output matrix.
         type(Diagonal) :: B
      end function construct_dense_to_diag
   end interface





   interface Bidiagonal
      pure module function initialize_bidiag(n, which) result(A)
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Upper- or lower-bidiagonal.
         character(len=1), optional, intent(in) :: which
         ! Output matrix.
         type(Bidiagonal) :: A
      end function initialize_bidiag

      pure module function construct_bidiag(dv, ev, which) result(A)
         ! Diagonal elements.
         real(wp), intent(in) :: dv(:), ev(:)
         ! Upper- or lower-bidiagonal.
         character(len=1), optional, intent(in) :: which
         ! Output matrix.
         type(Bidiagonal) :: A
      end function construct_bidiag

      pure module function construct_constant_bidiag(d, e, n, which) result(A)
         ! Diagonal elements.
         real(wp), intent(in) :: d, e
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Upper- or lower-bidiagonal.
         character(len=1), optional, intent(in) :: which
         ! Output matrix.
         type(Bidiagonal) :: A
      end function construct_constant_bidiag
   end interface




   interface Tridiagonal
      pure module function initialize_tridiag(n) result(A)
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Output matrix.
         type(Tridiagonal) :: A
      end function initialize_tridiag

      pure module function construct_tridiag(dl, d, du) result(A)
         ! Diagonals of the matrix.
         real(wp), intent(in) :: dl(:), d(:), du(:)
         ! Output matrix.
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
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Output matrix.
         type(SymTridiagonal) :: A
      end function initialize_symtridiag

      pure module function construct_symtridiag(dv, ev) result(A)
         ! Diagonal elements.
         real(wp), intent(in) :: dv(:), ev(:)
         ! Output matrix.
         type(SymTridiagonal) :: A
      end function construct_symtridiag

      pure module function construct_constant_symtridiag(d, e, n) result(A)
         ! Diagonal elements.
         real(wp), intent(in) :: d, e
         ! Dimension of the matrix.
         integer(int32), intent(in) :: n
         ! Output matrix.
         type(SymTridiagonal) :: A
      end function construct_constant_symtridiag
   end interface

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   interface matmul
      pure module function diag_spmv(A, x) result(y)
         ! Input matrix.
         type(Diagonal), intent(in) :: A
         ! Input vector.
         real(wp), intent(in) :: x(:)
         ! Output vector.
         real(wp) :: y(size(x))
      end function diag_spmv

      pure module function diag_multi_spmv(A, X) result(Y)
         ! Input matrix.
         type(Diagonal), intent(in) :: A
         ! Input vectors.
         real(wp), intent(in) :: X(:, :)
         ! Output vectors.
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function diag_multi_spmv

      pure module function bidiag_spmv(A, x) result(y)
         ! Input matrix.
         type(Bidiagonal), intent(in) :: A
         ! Input vector.
         real(wp), intent(in) :: x(:)
         ! Output vector.
         real(wp) :: y(size(x))
      end function bidiag_spmv

      pure module function bidiag_multi_spmv(A, X) result(Y)
         ! Input matrix.
         type(Bidiagonal), intent(in) :: A
         ! Input vectors.
         real(wp), intent(in) :: X(:, :)
         ! Output vectors.
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function bidiag_multi_spmv

      pure module function tridiag_spmv(A, x) result(y)
         ! Input matrix.
         type(Tridiagonal), intent(in) :: A
         ! Input vector.
         real(wp), intent(in) :: x(:)
         ! Output vector
         real(wp) :: y(size(x))
      end function tridiag_spmv

      pure module function tridiag_multi_spmv(A, X) result(Y)
         ! Input matrix.
         type(Tridiagonal), intent(in) :: A
         ! Input vectors.
         real(wp), intent(in) :: X(:, :)
         ! Output vectors.
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function tridiag_multi_spmv

      pure module function symtridiag_spmv(A, x) result(y)
         ! Input matrix.
         type(SymTridiagonal), intent(in) :: A
         ! Input vector.
         real(wp), intent(in) :: x(:)
         ! Output vector.
         real(wp) :: y(size(x))
      end function symtridiag_spmv

      pure module function symtridiag_multi_spmv(A, X) result(Y)
         ! Input matrix.
         type(SymTridiagonal), intent(in) :: A
         ! Input vectors.
         real(wp), intent(in) :: X(:, :)
         ! Output vectors.
         real(wp) :: Y(size(X, 1), size(X, 2))
      end function symtridiag_multi_spmv
   end interface

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   interface solve
      ! Diagonal matrix solve.
      pure module function diag_solve(A, b) result(x)
         type(Diagonal), intent(in) :: A
         real(wp), intent(in) :: b(:)
         real(wp) :: x(size(b))
      end function diag_solve

      pure module function diag_multi_solve(A, B) result(X)
         type(Diagonal), intent(in) :: A
         real(wp), intent(in) :: B(:, :)
         real(wp) :: X(size(B, 1), size(B, 2))
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
      pure module subroutine diag_eig(A, lambda, vectors)
         ! Input matrix.
         type(Diagonal), intent(in) :: A
         ! Eigenvalues.
         real(wp), intent(out) :: lambda(A%n)
         ! Eigenvectors.
         real(wp), intent(out) :: vectors(A%n, A%n)
      end subroutine diag_eig
   end interface

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
      pure module function diag_to_dense(A) result(B)
         ! Input Diagonal matrix.
         type(Diagonal), intent(in) :: A
         ! Output dense matrix.
         real(wp) :: B(A%n, A%n)
      end function diag_to_dense

      pure module function bidiag_to_dense(A) result(B)
         ! Input Bidiagonal matrix.
         type(Bidiagonal), intent(in) :: A
         ! Output dense matrix.
         real(wp) :: B(A%n, A%n)
      end function bidiag_to_dense

      pure module function tridiag_to_dense(A) result(B)
         ! Input tridiagonal matrix.
         type(Tridiagonal), intent(in) :: A
         ! Output dense matrix.
         real(wp) :: B(A%n, A%n)
      end function tridiag_to_dense

      pure module function symtridiag_to_dense(A) result(B)
         ! Input tridiagonal matrix.
         type(SymTridiagonal), intent(in) :: A
         ! Output dense matrix.
         real(wp) :: B(A%n, A%n)
      end function symtridiag_to_dense
   end interface





   interface transpose
      pure module function diag_transpose(A) result(B)
         ! Input diagonal matrix.
         type(Diagonal), intent(in) :: A
         ! Output diagonal matrix.
         type(Diagonal) :: B
      end function diag_transpose

      pure module function bidiag_transpose(A) result(B)
         ! Input bidiagonal matrix.
         type(Bidiagonal), intent(in) :: A
         ! Transposed matrix.
         type(Bidiagonal) :: B
      end function bidiag_transpose

      pure module function tridiag_transpose(A) result(B)
         ! Input tridiagonal matrix.
         type(Tridiagonal), intent(in) :: A
         ! Transposed matrix.
         type(Tridiagonal) :: B
      end function tridiag_transpose

      pure module function symtridiag_transpose(A) result(B)
         ! Input tridiagonal matrix.
         type(SymTridiagonal), intent(in) :: A
         ! Transposed matrix.
         type(SymTridiagonal) :: B
      end function symtridiag_transpose
   end interface

contains

end module SpecialMatrices_Tridiagonal
