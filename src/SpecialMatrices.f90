module SpecialMatrices
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
   public :: matmul
   public :: solve

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   public :: dense
   public :: shape

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, SpecialMatrices!"
   end subroutine say_hello
end module SpecialMatrices
