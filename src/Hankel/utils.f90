submodule(specialmatrices_hankel) hankel_utilities
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: i, j, m, n
   m = A%m; n = A%n; allocate (B(m, n), source=0.0_dp)
   do concurrent(j=1:n)
      B(:, j) = A%v(j:j + m - 1)
   end do
   end procedure dense_rdp

   module procedure transpose_rdp
   B = Hankel(A%v, A%n, A%m)
   end procedure transpose_rdp

   module procedure size_rdp
   if (present(dim)) then
      select case (dim)
      case (1)
         arr_size = A%m
      case (2)
         arr_size = A%n
      case default
         error stop "Matrix has only two dimensions."
      end select
   else
      arr_size = A%m*A%n
   end if
   end procedure size_rdp

   module procedure shape_rdp
   arr_shape = [A%m, A%n]
   end procedure shape_rdp

   module procedure scalar_multiplication_rdp
   B = Hankel(alpha*A%v, A%m, A%n)
   end procedure scalar_multiplication_rdp

   module procedure scalar_multiplication_bis_rdp
   B = Hankel(alpha*A%v, A%m, A%n)
   end procedure scalar_multiplication_bis_rdp
end submodule hankel_utilities
