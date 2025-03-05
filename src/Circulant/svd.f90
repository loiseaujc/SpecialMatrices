submodule (specialmatrices_circulant) circulant_svd
   use stdlib_sorting, only: sort, sort_index
   use stdlib_linalg, only: diag, hermitian
   implicit none(type, external)
contains
   module procedure svdvals_rdp
      !> Analytic expression for the singular values.
      s = abs(A%c_hat)
      !> Sort in decreasing order.
      call sort(s, reverse=.true.)
   end procedure

   module procedure svd_rdp
      integer(ilp) :: i, n
      integer(ilp), allocatable :: indices(:)
      complex(dp), allocatable :: lambda(:), left(:, :), right(:, :)
      complex(dp), allocatable :: pi(:, :), D(:, :)
      !> Matrix dimension.
      n = A%n
      !> Allocate variables.
      allocate(lambda(n), left(n, n), right(n, n))
      !> Eigendecomposition of the Circulant matrix.
      call eig(A, lambda, right=right, left=left) ; left = hermitian(left)
      !> Singular values.
      s = abs(A%c_hat)
      !> Right singular vectors.
      vt = discrete_hartley_transform(left)
      !> Left singular vectors.
      pi = matmul(transpose(left), left) ; D = diag(A%c_hat/s)
      u = discrete_hartley_transform(left)
      u = matmul(u, real(D%re - matmul(pi, D%im)))
      !> Sort the SVD.
      allocate(indices(n)) ; call sort_index(s, indices, reverse=.true.)
      u = u(:, indices) ; vt = vt(indices, :)
   end procedure

   pure function discrete_hartley_transform(F)  result(H)
      complex(dp), intent(in) :: F(:, :)
      !! Input matrix.
      real(dp), allocatable :: H(:, :)
      !! Hartley transform of F.
      H = F%re + F%im
   end function
end submodule
