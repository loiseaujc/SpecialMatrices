submodule(specialmatrices_diagonal) utilities
   use stdlib_linalg, only: diag
   implicit none(type, external)
contains
   module procedure dense_rdp
   B = diag(A%dv)
   end procedure

   module procedure transpose_rdp
   B = A
   end procedure

   module procedure size_rdp
   arr_size = A%n
   end procedure

   module procedure shape_rdp
   arr_shape = A%n
   end procedure

   module procedure scalar_multiplication_rdp
   B = Diagonal(alpha*A%dv)
   end procedure

   module procedure scalar_multiplication_bis_rdp
   B = Diagonal(alpha*A%dv)
   end procedure

end submodule
