submodule(specialmatrices_diagonal) diagonal_utilities
   use stdlib_linalg, only: diag
   implicit none(type, external)
contains
   module procedure dense_rdp
   B = diag(A%dv)
   end procedure dense_rdp

   module procedure transpose_rdp
   B = A
   end procedure transpose_rdp

   module procedure size_rdp
   arr_size = A%n
   end procedure size_rdp

   module procedure shape_rdp
   arr_shape = A%n
   end procedure shape_rdp

   module procedure scalar_multiplication_rdp
   B = Diagonal(alpha*A%dv)
   end procedure scalar_multiplication_rdp

   module procedure scalar_multiplication_bis_rdp
   B = Diagonal(alpha*A%dv)
   end procedure scalar_multiplication_bis_rdp

end submodule diagonal_utilities
