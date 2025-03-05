submodule(specialmatrices_circulant) circulant_utilities
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: j, n
   n = A%n ; allocate(B(n, n)) ; B = 0.0_dp
   do concurrent(j=1:n)
      B(:, j) = cshift(A%c, -j+1)
   enddo
   end procedure

   module procedure transpose_rdp
   real(dp), allocatable :: b_vec(:)
   integer(ilp) :: n
   n = A%n
   b_vec = A%c(n:1:-1) ; b_vec = cshift(b_vec, -1)
   B = Circulant(b_vec)
   end procedure

   module procedure size_rdp
   arr_size = A%n
   end procedure

   module procedure shape_rdp
   arr_shape = A%n
   end procedure

   module procedure scalar_multiplication_rdp
   B = Circulant(alpha*A%c)
   end procedure

   module procedure scalar_multiplication_bis_rdp
   B = Circulant(alpha*A%c)
   end procedure

end submodule
