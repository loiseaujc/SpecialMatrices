submodule(specialmatrices_hankel) hankel_linear_solver
    use specialmatrices_toeplitz, only: toeplitz_solve => solve

    implicit none(type, external)

contains
    module procedure solve_single_rhs
        integer(ilp) :: n
        type(Toeplitz) :: T
        real(dp), allocatable :: y(:)

        ! Matrix must be square.

        if (size(A, 1) /= size(A,2)) error stop "Matrix is not square."
        
        n = size (A, 2)
        
        if (size(b) /= n) error stop "Dimension of b is inconsistent with A"
        
        ! Hankel -> Toeplitz
        T = Toeplitz(A%v(n:), A%v(n:1:-1))

        ! Solve Ty = b
        y = toeplitz_solve(T, b)

        ! y -> x
        x = y(n:1:-1)

    end procedure solve_single_rhs

    module procedure solve_multi_rhs
        integer(ilp) :: n
        type(Toeplitz) :: T
        real(dp), allocatable :: Y(:, :)

        if (size(A, 1) /= size(A,2)) error stop "Matrix A is not square."
        
        n = size (A, 2)
        
        if (size(B, 1) /= n) error stop "Dimension of B is inconsistent with A"

        ! Hankel -> Toeplitz
        T = Toeplitz(A%v(n:), A%v(n:1:-1))
        
        ! Solve TY = B.
        Y = toeplitz_solve(T, B)

        ! Y -> X
        X = Y(n:1:-1, :)

    end procedure solve_multi_rhs

end submodule








