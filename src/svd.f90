submodule (specialmatrices_toeplitz) toeplitz_svd
   implicit none(type, external)
contains
   module procedure svdvals_rdp
      s = svdvals(dense(A))
   end procedure

   module procedure svd_rdp
      real(dp), allocatable :: Amat(:, :)
      Amat = dense(A) ; call svd(Amat, s, u, vt, overwrite_a=.true.)
   end procedure
end submodule
