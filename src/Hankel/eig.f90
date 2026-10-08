submodule(specialmatrices_hankel) hankel_eigendecomposition
    implicit none(type, external)
contains
    module procedure eigvalsh_rdp

        if (A%m /= A%n) error stop "Matrix A is not square"

        lambda = eigvalsh(dense(A))
    end procedure eigvalsh_rdp

    module procedure eigh_rdp
        real(dp), allocatable :: Amat(:, :)

        if (A%m /= A%n) error stop "Matrix A is not square"

        Amat = dense(A)
        call eigh(Amat, lambda, vectors=vectors, overwrite_a=.true.)
    end procedure eigh_rdp
end submodule hankel_eigendecomposition