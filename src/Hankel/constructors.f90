submodule(specialmatrices_hankel) hankel_constructors
   implicit none(type, external)
contains
   module procedure construct
   !> Sanity check.
   if (size(v) < m + n - 1) then
      error stop "Dimension of v is inconsistent with that of H."
   end if
   !> Dimension of the matrix.
   A%m = m; A%n = n
   !> Generating vector.
   A%v = v(1:m + n - 1)
   end procedure construct
end submodule hankel_constructors
