submodule(specialmatrices_toeplitz) toeplitz_utilities
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: i, j, m, n
   real(dp), allocatable :: t(:)
   m = A%m ; n = A%n ; allocate(t(-(n-1):m-1)) ; allocate(B(m, n))
   t(:-1) = A%vr(n:2:-1)
   t(0:) = A%vc
   do concurrent(i=1:m, j=1:n)
      B(i, j) = t(i-j)
   enddo
   end procedure

   module procedure transpose_rdp
   end procedure

   module procedure size_rdp
   if (present(dim)) then
      select case(dim)
         case (1)
            arr_size = A%m
         case (2)
            arr_size = A%n
      end select
   else
      arr_size = A%m * A%n
   endif
   end procedure

   module procedure shape_rdp
   arr_shape = [A%m, A%n]
   end procedure

   module procedure scalar_multiplication_rdp
   end procedure

   module procedure scalar_multiplication_bis_rdp
   end procedure

end submodule
