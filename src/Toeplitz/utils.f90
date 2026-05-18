submodule(specialmatrices_toeplitz) toeplitz_utilities
   implicit none(type, external)
contains
   module procedure dense_rdp
   integer(ilp) :: i, j, m, n
   real(dp), allocatable :: t(:)
   m = A%m; n = A%n; allocate (t(-(n - 1):m - 1)); allocate (B(m, n))
   t(:-1) = A%vr(n:2:-1)
   t(0:) = A%vc
   do concurrent(i=1:m, j=1:n)
      B(i, j) = t(i - j)
   end do
   end procedure dense_rdp

   module procedure transpose_rdp
   B = Toeplitz(A%vr, A%vc)
   end procedure transpose_rdp

   module procedure size_rdp
   if (present(dim)) then
      select case (dim)
      case (1)
         arr_size = A%m
      case (2)
         arr_size = A%n
      case default
         error stop "Matrix has only two dimensions."
      end select
   else
      arr_size = A%m*A%n
   end if
   end procedure size_rdp

   module procedure shape_rdp
   arr_shape = [A%m, A%n]
   end procedure shape_rdp

   module procedure scalar_multiplication_rdp
   B = Toeplitz(alpha*A%vc, alpha*A%vr)
   end procedure scalar_multiplication_rdp

   module procedure scalar_multiplication_bis_rdp
   B = Toeplitz(alpha*A%vc, alpha*A%vr)
   end procedure scalar_multiplication_bis_rdp

   module procedure Toeplitz2Circulant
   real(dp), allocatable :: c_vec(:)
   integer(ilp) :: m, n
   m = T%m; n = T%n; allocate (c_vec(m + n))
   c_vec(:m) = T%vc; c_vec(m + 1:) = cshift(T%vr(n:1:-1), -1)
   C = Circulant(c_vec)
   end procedure Toeplitz2Circulant
end submodule toeplitz_utilities
