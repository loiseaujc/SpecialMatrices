submodule(specialmatrices_symtridiagonal) symtridiagonal_matvecs
   implicit none(type, external)
contains
   module procedure spmv
   integer(ilp) :: i, n
   n = A%n; allocate (y, mold=x)
   y(1) = A%dv(1)*x(1) + A%ev(1)*x(2)
   do concurrent(i=2:n - 1)
      y(i) = A%ev(i - 1)*x(i - 1) + A%dv(i)*x(i) + A%ev(i)*x(i + 1)
   enddo
   y(n) = A%dv(n)*x(n) + A%ev(n - 1)*x(n - 1)
   end procedure

   module procedure spmvs
   integer(ilp) :: i, j, n, nvecs
   n = A%n; nvecs = size(x, 2); allocate (y, mold=x)
   y(1, :) = A%dv(1)*x(1, :) + A%ev(1)*x(2, :)
   do concurrent(i=2:n - 1, j=1:nvecs)
      y(i, j) = A%ev(i - 1)*x(i - 1, j) + A%dv(i)*x(i, j) + A%ev(i)*x(i + 1, j)
   enddo
   y(n, :) = A%dv(n)*x(n, :) + A%ev(n - 1)*x(n - 1, :)
   end procedure
end submodule
