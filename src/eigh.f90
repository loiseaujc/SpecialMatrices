submodule(specialmatrices_strang) strang_eigh
   use stdlib_constants, only: pi => pi_dp
   implicit none(type, external)
contains
   module procedure eigvalsh_rdp
   integer(ilp) :: k, n
   n = A%n; lambda = 2*[(1 - cos((pi*k)/(n + 1)), k=1, n)]
   end procedure

   module procedure eigh_rdp
   integer(ilp) :: i, j, n
   n = A%n
   lambda = eigvalsh(A)
   allocate (vectors(n, n))
   do concurrent(i=1:n, j=1:n)
      vectors(i, j) = sqrt(2.0_dp/(n + 1))*sin((i*j*pi)/(n + 1))
   end do
   end procedure
end submodule
