submodule(specialmatrices_poisson2D) poisson2D_solve
   use stdlib_linalg, only: eye, diag
   use stdlib_linalg_state, only: linalg_state_type, linalg_error_handling, LINALG_ERROR, &
                                  LINALG_INTERNAL_ERROR, LINALG_VALUE_ERROR, LINALG_SUCCESS
   use stdlib_io_npy, only: save_npy
   implicit none(type, external)

   character(len=*), parameter :: this = "poisson2D_linear_solver"
contains
   module procedure solve_single_rhs
      ! Local variables.
      real(dp), allocatable :: D(:), V(:, :)

      ! Compute eigendecomposition of the Poisson 2D.
      call eigh(A, D, V)

      ! Solution.
      x = matmul(transpose(V), b)
      x = matmul(diag(1.0_dp/D), x)
      x = matmul(V, x)
  end procedure

   module procedure solve_multi_rhs
      ! Local variables.
      real(dp), allocatable :: D(:), V(:, :)

      ! Compute eigendecomposition of the Poisson 2D.
      call eigh(A, D, V)

      ! Solution.
      x = matmul(transpose(V), b)
      x = matmul(diag(1.0_dp/D), x)
      x = matmul(V, x)
   end procedure

end submodule
