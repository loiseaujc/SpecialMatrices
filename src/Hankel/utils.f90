submodule(specialmatrices_hankel) hankel_utilities
   implicit none(type, external)
contains
   module procedure dense_rdp
      integer(ilp) :: i, j, m, n
      real(dp), allocatable :: t(:)
      m = A%m ; n = A%n ; allocate(t(-(n-1):m-1)) ; allocate(B(m, n))
      t(:-1) = A%vr(n:2:-1)
      t(0:) = A%vc
      do concurrent(i=1:m, j=1:n)
         B(i, j) = t(i-j)
      enddo
   end procedure dense_rdp

   module procedure transpose_rdp
      B = Hankel(A%vr, A%vc)
   end procedure transpose_rdp

   module procedure size_rdp
      if (present(dim)) then
         select case(dim)
            case (1)
               arr_size = A%m
            case (2)
               arr_size = A%n
            case default
               error stop "Matrix has only two dimensions."
         end select
      else
         arr_size = A%m * A%n
      endif
   end procedure size_rdp

   module procedure shape_rdp
      arr_shape = [A%m, A%n]
   end procedure shape_rdp

   module procedure scalar_multiplication_rdp
      B = Hankel(alpha*A%vc, alpha*A%vr)
   end procedure scalar_multiplication_rdp

   module procedure scalar_multiplication_bis_rdp
      B = Hankel(alpha*A%vc, alpha*A%vr)
   end procedure scalar_multiplication_bis_rdp

   module procedure Hankel2Toeplitz
      integer(ilp) :: m, n
      real(dp), allocatable :: vc(:), vr(:)
      m = H%m ; n = H%n
      vc = H%vr(2:2+m) ; vr = H%vc  ! Incorrect
      T = Toeplitz(vc, vr)
   end procedure Hankel2Toeplitz
end submodule hankel_utilities
