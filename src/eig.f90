submodule(specialmatrices_bidiagonal) bidiagonal_eigenvalue_decomposition
   use stdlib_linalg, only: stdlib_eig => eig, stdlib_eigvals => eigvals
   implicit none(type, external)
contains
   module procedure eigvals_rdp
   lambda = stdlib_eigvals(dense(A))
   end procedure

   module procedure eig_rdp
   real(dp), allocatable :: Amat(:, :)
   Amat = dense(A)
   call stdlib_eig(Amat, lambda, right=right, left=left, overwrite_a=.true.)
   end procedure
end submodule
