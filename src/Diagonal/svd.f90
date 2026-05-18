submodule(specialmatrices_diagonal) diagonal_singular_value_decomposition
   use stdlib_sorting, only: sort, sort_index
   use stdlib_linalg, only: eye
   implicit none(type, external)
contains
   module procedure svdvals_rdp
   s = abs(A%dv); call sort(s, reverse=.true.)
   end procedure svdvals_rdp

   module procedure svd_rdp
   integer(ilp) :: i, index(A%n)
   ! Sorted singular values.
   s = abs(A%dv); call sort_index(s, index, reverse=.true.)
   ! Left singular vectors.
   u = eye(A%n); u = u(:, index)
   ! Right singular vectors.
   vt = eye(A%n)
   do concurrent(i=1:A%n, A%dv(i) < 0.0_dp)
      vt(i, i) = -1.0_dp
   end do
   vt = vt(index, :)
   end procedure svd_rdp
end submodule diagonal_singular_value_decomposition
