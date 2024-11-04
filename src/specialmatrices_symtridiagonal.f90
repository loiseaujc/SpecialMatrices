submodule(SpecialMatrices_Tridiagonal) SymmetricTridiagonal
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: gtsv, ptsv
   implicit none(type, external)

contains

   !--------------------------------
   !-----     Constructors     -----
   !--------------------------------

   module procedure initialize_symtridiag
   A%n = n; allocate (A%dv(n), A%ev(n - 1))
   A%dv = 0.0_dp; A%ev = 0.0_dp; A%isposdef = .false.
   end procedure initialize_symtridiag

   module procedure construct_symtridiag
   A%n = size(dv); A%dv = dv; A%ev = ev; A%isposdef = optval(isposdef, .false.)
   end procedure construct_symtridiag

   module procedure construct_constant_symtridiag
   integer(ilp) :: i
   A%n = n; A%dv = [(d, i=1, n)]; A%ev = [(e, i=1, n - 1)]; A%isposdef = optval(isposdef, .false.)
   end procedure construct_constant_symtridiag

   !------------------------------------------------------------
   !-----     Matrix-vector and Matrix-matrix products     -----
   !------------------------------------------------------------

   module procedure symtridiag_trace
   tr = sum(A%dv)
   end procedure symtridiag_trace

   module procedure symtridiag_det
   real(dp) :: f_0, f_1
   integer(ilp) :: i
   f_0 = 1.0_dp; f_1 = 0.0_dp
   ! First iteration.
   determinant = A%dv(1)*f_0; f_1 = f_0; f_0 = determinant
   ! Continuants
   do i=2, A%n
      ! Recurrence relation.
      determinant = A%dv(i)*f_0 - A%ev(i-1)**2 * f_1
      ! Store previous values.
      f_1 = f_0; f_0 = determinant
   enddo
   end procedure symtridiag_det

   module procedure symtridiag_spmv
   integer(ilp) :: i, n
   n = size(x); allocate (y, mold=x)
   y(1) = A%dv(1)*x(1) + A%ev(1)*x(2)
   do concurrent(i=2:n - 1)
      y(i) = A%ev(i - 1)*x(i - 1) + A%dv(i)*x(i) + A%ev(i)*x(i + 1)
   end do
   y(n) = A%dv(n)*x(n) + A%ev(n - 1)*x(n - 1)
   end procedure symtridiag_spmv

   module procedure symtridiag_multi_spmv
   integer(ilp) :: i, j, n
   n = size(x, 1); allocate (y, mold=x)
   y(1, :) = A%dv(1)*x(1, :) + A%ev(1)*x(2, :)
   do concurrent(i=2:n - 1, j=1:size(x, 2))
      y(i, j) = A%ev(i - 1)*x(i - 1, j) + A%dv(i)*x(i, j) + A%ev(i)*x(i + 1, j)
   end do
   y(n, :) = A%dv(n)*x(n, :) + A%ev(n - 1)*x(n - 1, :)
   end procedure symtridiag_multi_spmv

   module procedure symtridiag_spmv_ip
   integer(ilp) :: i, n
   n = size(x)
   y(1) = A%dv(1)*x(1) + A%ev(1)*x(2)
   do concurrent(i=2:n - 1)
      y(i) = A%ev(i - 1)*x(i - 1) + A%dv(i)*x(i) + A%ev(i)*x(i + 1)
   end do
   y(n) = A%dv(n)*x(n) + A%ev(n - 1)*x(n - 1)
   return
   end procedure symtridiag_spmv_ip

   module procedure symtridiag_multi_spmv_ip
   integer(ilp) :: i, j, n
   n = size(x, 1)
   y(1, :) = A%dv(1)*x(1, :) + A%ev(1)*x(2, :)
   do concurrent(i=2:n - 1, j=1:size(x, 2))
      y(i, j) = A%ev(i - 1)*x(i - 1, j) + A%dv(i)*x(i, j) + A%ev(i)*x(i + 1, j)
   end do
   y(n, :) = A%dv(n)*x(n, :) + A%ev(n - 1)*x(n - 1, :)
   return
   end procedure symtridiag_multi_spmv_ip

   module procedure symtridiag_eigvalsh
      call symtridiag_eigh_driver(A, lambda)
   end procedure symtridiag_eigvalsh

   module procedure symtridiag_eigh
      call symtridiag_eigh_driver(A, lambda, vectors)
   end procedure symtridiag_eigh

   subroutine symtridiag_eigh_driver(A, lambda, vectors)
      type(SymTridiagonal), intent(in) :: A
      real(dp), allocatable, intent(out) :: lambda(:)
      real(dp), allocatable, optional, intent(out) :: vectors(:, :)

      ! Local variables.
      return
   end subroutine symtridiag_eigh_driver

   !----------------------------------
   !-----     Linear Algebra     -----
   !----------------------------------

   module procedure symtridiag_solve
   integer(ilp) :: i, n, nrhs, info
   real(dp) :: dl(A%n - 1), d(A%n), du(A%n - 1), b_(A%n, 1)
   ! Initialize arrays.
   n = A%n; dl = A%ev; d = A%dv; du = A%ev; b_(:, 1) = b; nrhs = 1
   ! Solve the system.
   if (A%isposdef) then
      call ptsv(n, nrhs, d, du, b_, n, info)
   else
      call gtsv(n, nrhs, dl, d, du, b_, n, info)
   end if
   ! Return the solution.
   x = b_(:, 1)
   end procedure symtridiag_solve

   module procedure symtridiag_multi_solve
   integer(ilp) :: i, n, nrhs, info
   real(dp) :: dl(A%n - 1), d(A%n), du(A%n - 1)
   ! Initialize arrays.
   n = A%n; dl = A%ev; d = A%dv; du = A%ev; nrhs = size(B, 2); X = B
   ! Solve the systems.
   if (A%isposdef) then
      call ptsv(n, nrhs, d, du, X, n, info)
   else
      call gtsv(n, nrhs, dl, d, du, X, n, info)
   end if
   end procedure symtridiag_multi_solve

   !-------------------------------------
   !-----     Utility procedures     -----
   !-------------------------------------

   module procedure symtridiag_to_dense
   integer(ilp) :: i, n
   n = A%n; B = 0.0_dp
   B(1, 1) = A%dv(1); B(1, 2) = A%ev(1)
   do concurrent(i=2:n - 1)
      B(i, i - 1) = A%ev(i - 1)
      B(i, i) = A%dv(i)
      B(i, i + 1) = A%ev(i)
   end do
   B(n, n - 1) = A%ev(n - 1); B(n, n) = A%dv(n)
   end procedure symtridiag_to_dense

   module procedure symtridiag_transpose
   B = A
   end procedure symtridiag_transpose

   module procedure symtridiag_shape
   shape(2) = A%n
   end procedure symtridiag_shape
   
   module procedure symtridiag_size
   ! Utility function to get the size of a `SymTridiagonal` matrix along a given dimension.
   arr_size = A%n
   end procedure symtridiag_size

   module procedure symtridiag_scalar_multiplication
   B = A; B%dv = alpha*A%dv; B%ev = alpha*A%ev
   if (alpha <= 0.0_dp) B%isposdef = .false.
   end procedure symtridiag_scalar_multiplication

   module procedure symtridiag_scalar_multiplication_bis
   B = symtridiag_scalar_multiplication(alpha, A)
   end procedure symtridiag_scalar_multiplication_bis

end submodule SymmetricTridiagonal
