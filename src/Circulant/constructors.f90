submodule(specialmatrices_circulant) circulant_constructors
   implicit none(type, external)
contains
   module procedure initialize
   A%n = n; allocate(A%c(n)) ; A%c = 0.0_dp
   !> Initialize the workspace for the FFT and its inverse.
   allocate(A%wsave(2*n+15)) ; call init_fft(n, A%wsave)
   end procedure

   module procedure construct
   integer(ilp) :: n
   !> Initialize the standard matrix data.
   n = size(c) ; A%n = n; A%c = c
   !> Initialize the workspace for the FFT and its inverse.
   allocate(A%wsave(2*n+15)) ; call init_fft(n, A%wsave)
   end procedure
end submodule
