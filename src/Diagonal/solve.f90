submodule(specialmatrices_diagonal) diagonal_linear_solver
   implicit none(type, external)
contains
   module procedure solve_single_rhs
   ! Solve \(Ax = b\).
   integer(ilp) :: i
   allocate (x, mold=b)
   do concurrent(i=1:A%n)
      x(i) = b(i)/A%dv(i)
   end do
   end procedure

   module procedure solve_multi_rhs
   ! Solve \(AX = B\).
   integer(ilp) :: i, j
   real(dp), allocatable :: inv_dv(:)
   allocate (x, mold=b); inv_dv = 1.0_dp/A%dv
   do concurrent(i=1:A%n, j=1:size(b, 2))
      x(i, j) = b(i, j)*inv_dv(i)
   end do
   end procedure
end submodule
