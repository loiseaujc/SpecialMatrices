module SpecialMatrices
   use SpecialMatrices_Tridiagonal
   implicit none
   private

   !--------------------------------
   !-----     Matrix types     -----
   !--------------------------------

   public :: Tridiagonal

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   public :: matmul
   public :: solve

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   public :: transpose
   public :: dense

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, SpecialMatrices!"
   end subroutine say_hello
end module SpecialMatrices
