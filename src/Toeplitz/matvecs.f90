submodule(specialmatrices_toeplitz) toeplitz_matvecs
contains
   module procedure spmv
   integer(ilp) :: m, n
   real(dp), allocatable :: x_circ(:), y_circ(:)
   !> Dimension of the matrix.
   m = A%m ; n = A%n
   !> Allocate variables.
   allocate(x_circ(m+n)) ; x_circ = 0.0_dp ; x_circ(:n) = x
   !> Toeplitz spmv via Circulant embedding.
   y_circ = matmul(Circulant(A), x_circ) ; y = y_circ(:m)
   end procedure

   module procedure spmvs
   integer(ilp) :: m, n
   real(dp), allocatable :: x_circ(:, :), y_circ(:, :)
   !> Dimension of the matrix.
   m = A%m ; n = A%n
   !> Allocate variables.
   allocate(x_circ(m+n, size(x, 2))) ; x_circ = 0.0_dp ; x_circ(:n, :) = x
   !> Toeplitz spmv via Circulant embedding.
   y_circ = matmul(Circulant(A), x_circ) ; y = y_circ(:m, :)
   end procedure
end submodule
