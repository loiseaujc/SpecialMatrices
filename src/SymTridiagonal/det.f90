submodule(specialmatrices_symtridiagonal) symtridiagonal_determinant
   implicit none(type, external)
contains
   module procedure det_rdp
   real(dp) :: f_0, f_1
   integer(ilp) :: i
   f_0 = 1.0_dp; f_1 = 0.0_dp
   ! First iteration.
   d = A%dv(1)*f_0; f_1 = f_0; f_0 = d
   ! Continuants
   do i = 2, A%n
      ! Reccurence relation.
      d = A%dv(i)*f_0 - A%ev(i - 1)**2*f_1
      ! Store previous values.
      f_1 = f_0; f_0 = d
   end do
   end procedure
end submodule
