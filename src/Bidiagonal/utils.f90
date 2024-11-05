submodule(specialmatrices_bidiagonal) bidiagonal_utils
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: i, n
   n = A%n; allocate (B(n, n)); B = 0.0_dp
   select case(A%which)
      case("L")
         B(1, 1) = A%dv(1)
         do concurrent(i=2:n)
            B(i, i) = A%dv(i)
            B(i, i - 1) = A%ev(i - 1)
         enddo
      case("U")
      do concurrent(i=1:n - 1)
         B(i, i) = A%dv(i)
         B(i, i + 1) = A%ev(i)
      enddo
      B(n, n) = A%dv(n)
   end select
   end procedure

   module procedure transpose_rdp
   B = A
   if (A%which == "L") then
      B%which = "U"
   else
      B%which = "L"
   endif
   end procedure

   module procedure shape_rdp
   arr_shape = A%n
   end procedure

   module procedure size_rdp
   arr_size = A%n
   end procedure

   module procedure scalar_multiplication_rdp
   B = Bidiagonal(alpha*A%dv, alpha*A%ev, A%which)
   end procedure

   module procedure scalar_multiplication_bis_rdp
   B = Bidiagonal(alpha*A%dv, alpha*A%ev, A%which)
   end procedure

end submodule
