submodule (specialmatrices_toeplitz) toeplitz_eigendecomposition
   implicit none(type, external)
contains
   module procedure eigvals_rdp
      lambda = eigvals(dense(A))
   end procedure

   module procedure eig_rdp
      real(dp), allocatable :: Amat(:, :)
      Amat = dense(A) ; call eig(Amat, lambda, right=right, left=left, overwrite_a=.true.)
   end procedure
end submodule
