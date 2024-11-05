submodule(specialmatrices_diagonal) diagonal_constructors
   use stdlib_linalg, only: diag
   implicit none(type, external)
contains
   module procedure initialize
   ! Utility function to construct a `Diagonal` matrix filled with zeros.
   A%n = n; allocate (A%dv(n)); A%dv = 0.0_dp
   end procedure

   module procedure construct
   ! Utility function to construct a `Diagonal` matrix from a rank-1 array.
   A%n = size(dv); A%dv = dv
   end procedure

   module procedure construct_constant
   ! Utility function to construct a `Diagonal` matrix with constant value.
   integer(ilp) :: i
   A%n = n; A%dv = [(d, i=1, n)]
   end procedure

   module procedure dense_to_diag
   ! Utility function to construct a `Diagonal` matrix from a rank-2 array.
   B = Diagonal(diag(A))
   end procedure
end submodule
