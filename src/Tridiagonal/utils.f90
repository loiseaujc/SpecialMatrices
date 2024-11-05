submodule(specialmatrices_tridiagonal) tridiagonal_utils
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: i, n
   n = A%n; allocate (B(n, n)); B = 0.0_dp
   B(1, 1) = A%dv(1); B(1, 2) = A%du(1)
   do concurrent(i=2:n - 1)
      B(i, i - 1) = A%dl(i - 1)
      B(i, i) = A%dv(i)
      B(i, i + 1) = A%du(i)
   end do
   B(n, n - 1) = A%dl(n - 1); B(n, n) = A%dv(n)
   end procedure

   module procedure transpose_rdp
   B = A
   end procedure

   module procedure shape_rdp
   arr_shape = A%n
   end procedure

   module procedure size_rdp
   arr_size = A%n
   end procedure

   module procedure scalar_multiplication_rdp
   B = Tridiagonal(alpha*A%dl, alpha*A%dv, alpha*A%du)
   end procedure

   module procedure scalar_multiplication_bis_rdp
   B = Tridiagonal(alpha*A%dl, alpha*A%dv, alpha*A%du)
   end procedure

end submodule
