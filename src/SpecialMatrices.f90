module SpecialMatrices
   use specialmatrices_diagonal
   use specialmatrices_symtridiagonal
   use SpecialMatrices_Tridiagonal
   implicit none
   private

   !--------------------------------
   !-----     Matrix types     -----
   !--------------------------------

   public :: Diagonal
   public :: Bidiagonal
   public :: Tridiagonal
   public :: SymTridiagonal

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   public :: transpose
   public :: det
   public :: trace
   public :: inv
   public :: matmul
   public :: solve
   public :: svd, svdvals
   public :: eigh, eigvalsh

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, SpecialMatrices!"
   end subroutine say_hello
end module SpecialMatrices
