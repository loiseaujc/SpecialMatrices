submodule(specialmatrices_diagonal) diagonal_hermitian_eigenvalue_decomposition
   use stdlib_sorting, only: sort, sort_index
   use stdlib_linalg, only: eye
   implicit none(type, external)
contains
   module procedure eigvalsh_rdp
   lambda = A%dv; call sort(lambda)
   end procedure

   module procedure eigh_rdp
   integer(ilp) :: index(A%n)
   ! Eigenvalues.
   lambda = A%dv; call sort_index(lambda, index)
   ! Eigenvectors.
   if (present(vectors)) then
      vectors = eye(A%n); vectors = vectors(:, index)
   end if
   end procedure
end submodule
