submodule(specialmatrices_circulant) circulant_linear_solver
   implicit none(type, external)
contains
   module procedure solve_single_rhs
      x = real(ifft(fft(cmplx(b, kind=dp), A%n) / A%c_hat, A%n)) / A%n
   end procedure

   module procedure solve_multi_rhs
      integer(ilp) :: i
      x = b
      do concurrent(i=1:size(b, 2))
         x(:, i) = solve(A, b(:, i))
      enddo
   end procedure
end submodule
