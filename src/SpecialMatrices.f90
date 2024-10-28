module SpecialMatrices
   use SpecialMatrices_Tridiagonal
   implicit none
   private

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, SpecialMatrices!"
   end subroutine say_hello
end module SpecialMatrices
