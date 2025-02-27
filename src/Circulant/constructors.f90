submodule(specialmatrices_circulant) circulant_constructors
   implicit none(type, external)
contains
   module procedure initialize
   A%n = n; allocate(A%c(n)) ; A%c = 0.0_dp
   !> Initialize the workspace for the FFT and its inverse.
   A%c_hat = 0.0_dp
   end procedure

   module procedure construct
   integer(ilp) :: n
   !> Initialize the standard matrix data.
   n = size(c) ; A%n = n; A%c = c
   !> Fourier Transform of the generating vector.
   A%c_hat = fft(cmplx(c, kind=dp), n)
   end procedure
end submodule
