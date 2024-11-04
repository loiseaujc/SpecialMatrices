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
   public :: det
   public :: trace
   public :: inv
   public :: matmul, spmv_ip
   public :: solve, solve_ip
   public :: svd, svdvals

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   public :: dense
   public :: shape
   public :: size

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, SpecialMatrices!"
   end subroutine say_hello
end module SpecialMatrices
