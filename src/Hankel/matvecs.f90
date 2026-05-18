submodule(specialmatrices_hankel) hankel_matvecs
   implicit none(type, external)
contains
   module procedure spmv
   integer(ilp) :: m, n
   m = A%m; n = A%n
   y = matmul(Toeplitz(A%v(n:), A%v(n:1:-1)), x(n:1:-1))
   end procedure spmv

   module procedure spmvs
   integer(ilp) :: m, n
   m = A%m; n = A%n
   y = matmul(Toeplitz(A%v(n:), A%v(n:1:-1)), x(n:1:-1, :))
   end procedure spmvs
end submodule hankel_matvecs
