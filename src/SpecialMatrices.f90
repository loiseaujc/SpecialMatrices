module SpecialMatrices
   use specialmatrices_diagonal
   use specialmatrices_bidiagonal
   use specialmatrices_tridiagonal
   use specialmatrices_symtridiagonal
   use specialmatrices_strang
   use specialmatrices_poisson2D
   use specialmatrices_circulant
   use specialmatrices_toeplitz
   use specialmatrices_hankel
   implicit none(type, external)
   private

   !--------------------------------
   !-----     Matrix types     -----
   !--------------------------------

   public :: Diagonal
   public :: Bidiagonal
   public :: Tridiagonal
   public :: SymTridiagonal
   public :: Strang
   public :: Poisson2D
   public :: Circulant
   public :: Toeplitz
   public :: Hankel

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
   public :: eig, eigvals

   !-------------------------------------
   !-----     Utility functions     -----
   !-------------------------------------

   public :: dense
   public :: shape
   public :: size
   public :: operator(*)

end module SpecialMatrices
