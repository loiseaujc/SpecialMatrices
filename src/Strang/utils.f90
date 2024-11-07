submodule(specialmatrices_strang) strang_utils
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: i, n
   n = A%n
   allocate (B(n, n)); B = 0.0_dp
   B(1, 1) = 2; B(1, 2) = -1
   do concurrent(i=2:A%n - 1)
      B(i, i - 1) = -1; B(i, i) = 2; B(i, i + 1) = -1
   end do
   B(n, n - 1) = -1; B(n, n) = 2
   end procedure

   module procedure shape_rdp
   arr_shape = A%n
   end procedure

   module procedure size_rdp
   arr_size = A%n
   end procedure
end submodule
