submodule(specialmatrices_circulant) circulant_constructors
   implicit none(type, external)
contains
   module procedure construct
   integer(ilp) :: n
   !> Initialize the standard matrix data.
   n = size(c) ; A%n = n; A%c = c
   !> Fourier Transform of the generating vector.
   A%c_hat = fft(cmplx(c, kind=dp), n)
   end procedure
end submodule
