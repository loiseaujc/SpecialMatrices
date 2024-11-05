submodule(specialmatrices_diagonal) diagonal_matvecs
   implicit none(type, external)
contains
   module procedure spmv
   ! Utility function to compute the matrix-vector product.
   integer(ilp) :: i
   allocate (y, mold=x)
   do concurrent(i=1:size(x))
      y(i) = A%dv(i)*x(i)
   end do
   end procedure

   module procedure spmvs
   ! Utility function to compute multiple matrix-vector products.
   integer(ilp) :: i, j
   allocate (y, mold=x)
   do concurrent(i=1:size(x, 1), j=1:size(x, 2))
      y(i, j) = A%dv(i)*x(i, j)
   end do
   end procedure
end submodule
