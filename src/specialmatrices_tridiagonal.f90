module SpecialMatrices_Tridiagonal
   use stdlib_kinds, only: wp => dp, int32
   implicit none
   private
   
   ! --> Constructors.
   public :: Tridiagonal
   public :: SymTridiagonal

   ! --> Linear Algebra.
   public :: matmul
   public :: solve
   public :: transpose

   ! --> Utility functions.
   public :: dense

   !-------------------------------------------------------
   !-----     Base types for tridiagonal matrices     -----
   !-------------------------------------------------------

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

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   interface dense
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
