submodule(specialmatrices_strang) strang_matvecs
   implicit none(type, external)
contains
   module procedure spmv
   integer(ilp) :: i, n
   n = A%n; allocate (y, mold=x)
   y(1) = 2*x(1) - x(2)
   do concurrent(i=2:n - 1)
      y(i) = -x(i - 1) + 2*x(i) - x(i + 1)
   end do
   y(n) = 2*x(n) - x(n - 1)
   end procedure

   module procedure spmvs
   integer(ilp) :: i
   allocate (y, mold=x)
   do concurrent(i=1:size(x, 2))
      y(:, i) = spmv(A, x(:, i))
   end do
   end procedure
end submodule
