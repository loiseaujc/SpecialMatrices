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
   end procedure

   module procedure size_rdp
   end procedure

   module procedure shape_rdp
   end procedure

   module procedure scalar_multiplication_rdp
   end procedure

   module procedure scalar_multiplication_bis_rdp
   end procedure

end submodule
