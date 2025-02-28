submodule(specialmatrices_toeplitz) toeplitz_matvecs
   use specialmatrices_circulant
contains
   module procedure spmv
   integer(ilp) :: m, n
   type(Circulant) :: C
   real(dp), allocatable :: c_vec(:)
   !> Dimension of the matrix.
   m = A%m ; n = A%n
   !> Circulant vector.
   allocate(c_vec(m+n))
   c_vec(:m) = A%vc ; c_vec(m+1:) = cshift(A%vr(n:1:-1), -1)
   !> Circulant matrix.
   C = Circulant(c_vec)
   !> Matrix-vector product.
   block
   real(dp), allocatable :: x_circ(:), y_circ(:)
   allocate(x_circ(m+n)) ; x_circ = 0.0_dp ; x_circ(:n) = x
   y_circ = matmul(C, x_circ) ; y = y_circ(:m)
   end block
   end procedure

   module procedure spmvs
   integer(ilp) :: m, n
   type(Circulant) :: C
   real(dp), allocatable :: c_vec(:)
   !> Dimension of the matrix.
   m = A%m ; n = A%n
   !> Circulant vector.
   allocate(c_vec(m+n))
   c_vec(:m) = A%vc ; c_vec(m+1:) = cshift(A%vr(n:1:-1), -1)
   !> Circulant matrix.
   C = Circulant(c_vec)
   !> Matrix-vector product.
   block
   real(dp), allocatable :: x_circ(:, :), y_circ(:, :)
   allocate(x_circ(m+n, size(x, 2))) ; x_circ = 0.0_dp ; x_circ(:n, :) = x
   y_circ = matmul(C, x_circ) ; y = y_circ(:m, :)
   end block
   end procedure
end submodule
