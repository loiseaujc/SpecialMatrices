submodule(specialmatrices_bidiagonal) bidiagonal_matvecs
   implicit none(type, external)
contains
   module procedure spmv
   integer(ilp) :: i, n
   n = A%n; allocate (y, mold=x)
   select case(A%which)
      case("L")
         y(1) = A%dv(1)*x(1)
         do concurrent(i=2:n)
            y(i) = A%ev(i - 1)*x(i - 1) + A%dv(i)*x(i)
         enddo
      case("U")
      do concurrent(i=1:n - 1)
         y(i) = A%dv(i)*x(i) + A%ev(i)*x(i + 1)
      enddo
      y(n) = A%dv(n)*x(n)
   end select
   end procedure

   module procedure spmvs
   integer(ilp) :: i, j, n, nvecs
   n = A%n; nvecs = size(x, 2); allocate (y, mold=x)
   do concurrent(i=1:nvecs)
      y(:, i) = spmv(A, x(:, i))
   enddo
   end procedure
end submodule
