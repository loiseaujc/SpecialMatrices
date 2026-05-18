submodule(specialmatrices_circulant) circulant_eigenvalue_decomposition
   use stdlib_linalg, only: hermitian, eye
   implicit none(type, external)
contains
   module procedure eigvals_rdp
   lambda = A%c_hat
   end procedure eigvals_rdp

   module procedure eig_rdp
   real(dp), allocatable :: Amat(:, :)
   integer(ilp) :: i, n
   n = A%n
   lambda = A%c_hat
   if (present(right)) then
      right = eye(n, mold=1.0_dp)
      do concurrent(i=1:n)
         right(:, i) = fft(right(:, i), n) / sqrt(1.0_dp*n)
      enddo
      right = hermitian(right)
   endif

   if (present(left)) then
      left = eye(n, mold=1.0_dp)
      do concurrent(i=1:n)
         left(:, i) = ifft(left(:, i), n) / sqrt(1.0_dp*n)
      enddo
   endif
   end procedure eig_rdp
end submodule circulant_eigenvalue_decomposition
