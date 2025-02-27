submodule(specialmatrices_circulant) circulant_matvecs
   implicit none(type, external)
contains
   module procedure spmv
      integer(ilp) :: i
      y = x
      y = real(ifft(fft(cmplx(y, kind=dp), A%n) * A%c_hat, A%n)) / A%n
   end procedure

   module procedure spmvs
      integer(ilp) :: i
      y = x
      do concurrent(i=1:size(x, 2))
      y(:, i) = matmul(A, x(:, i))
      enddo
   end procedure
end submodule
