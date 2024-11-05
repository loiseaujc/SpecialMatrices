submodule(specialmatrices_tridiagonal) tridiagonal_singular_value_decomposition
   use stdlib_linalg, only: stdlib_svd => svd, stdlib_svdvals => svdvals
   implicit none(type, external)
contains
   module procedure svdvals_rdp
   s = stdlib_svdvals(dense(A))
   end procedure

   module procedure svd_rdp
   real(dp), allocatable :: Amat(:, :)
   Amat = dense(A)
   call stdlib_svd(Amat, s, u, vt, overwrite_a=.true.)
   end procedure
end submodule
