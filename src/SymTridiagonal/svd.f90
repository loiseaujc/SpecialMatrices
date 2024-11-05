submodule(specialmatrices_symtridiagonal) symtridiagonal_singular_value_decomposition
   implicit none(type, external)
contains
   module procedure svdvals_rdp
   ! Get singular values from the eigendecomposition.
   s = abs(eigvalsh(A))
   end procedure

   module procedure svd_rdp
   integer(ilp) :: i, n
   n = A%n
   ! Compute eigendecomposition.
   call eigh(A, s, u)
   ! Left singular vectors.
   if (present(vt)) then
      vt = u
      do concurrent(i=1:n, s(i) < 0.0_dp)
         vt(:, i) = -1.0_dp*vt(:, i)
      end do
      vt = transpose(vt)
   end if
   ! Singular values.
   s = abs(s)
   end procedure
end submodule
