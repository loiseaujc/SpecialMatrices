module test_toeplitz
   ! Fortran standard library.
   use stdlib_math, only: is_close, all_close
   use stdlib_sorting, only: sort_index
   use stdlib_linalg_constants, only: dp, ilp
   use stdlib_linalg, only: diag, det, trace, inv, solve, svdvals, eigvals, hermitian
   ! Testdrive.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   ! SpecialMatrices
   use SpecialMatrices
   implicit none
   private

   integer, parameter :: m = 128, n = 64

   public :: collect_toeplitz_testsuite
contains

   !--------------------------------------
   !-----     toeplitz MATRICES     -----
   !--------------------------------------

   subroutine collect_toeplitz_testsuite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  ! new_unittest("toeplitz scalar multiplication", test_scalar_multiplication), &
                  ! new_unittest("toeplitz trace", test_trace), &
                  ! new_unittest("toeplitz determinant", test_det), &
                  ! new_unittest("toeplitz inverse", test_inv), &
                  new_unittest("toeplitz matmul", test_matmul) &
                  ! new_unittest("toeplitz linear solver", test_solve), &
                  ! new_unittest("toeplitz svdvals", test_svdvals), &
                  ! new_unittest("toeplitz svd", test_svd), &
                  ! new_unittest("toeplitz eigvals", test_eigvals), &
                  ! new_unittest("toeplitz eig", test_eig) &
                  ]
      return
   end subroutine collect_toeplitz_testsuite

   ! subroutine test_scalar_multiplication(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A, B
   !    real(dp), allocatable :: dv(:)
   !    real(dp) :: alpha
   !
   !    ! Initialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !    call random_number(alpha)
   !
   !    ! Scalar-matrix multiplication.
   !    B = alpha*A
   !    ! Check error.
   !    call check(error, all_close(alpha*dense(A), dense(B)), &
   !               "alpha*toeplitz failed.")
   !    if (allocated(error)) return
   !
   !    ! Matrix-scalar multipliation.
   !    B = A*alpha
   !    ! Check error.
   !    call check(error, all_close(alpha*dense(A), dense(B)), &
   !               "toeplitz*alpha failed.")
   !    return
   ! end subroutine test_scalar_multiplication
   !
   ! subroutine test_trace(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:)
   !
   !    ! Initialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !
   !    ! Compare against stdlib_linalg implementation.
   !    call check(error, is_close(trace(A), trace(dense(A))), &
   !               "toeplitz trace failed.")
   !    return
   ! end subroutine test_trace
   !
   ! subroutine test_det(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:)
   !
   !    ! Initialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !
   !    ! Compare against stdlib_linalg implementation.
   !    call check(error, is_close(det(A), det(dense(A))), &
   !               "toeplitz det failed.")
   !    return
   ! end subroutine test_det
   !
   ! subroutine test_inv(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:)
   !
   !    ! Intialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !
   !    ! Compare against stdlib_linalg implementation.
   !    call check(error, all_close(inv(dense(A)), dense(inv(A))), &
   !               "toeplitz inv failed.")
   !    return
   ! end subroutine test_inv

   subroutine test_matmul(error)
      type(error_type), allocatable, intent(out) :: error
      type(Toeplitz) :: A
      real(dp), allocatable :: vc(:), vr(:)

      ! Initialize matrix.
      allocate(vc(m)) ; call random_number(vc)
      allocate(vr(n)) ; call random_number(vr)
      A = Toeplitz(vc, vr)

      ! Matrix-vector product.
      block
         real(dp), allocatable :: x(:), y(:), y_dense(:)
         allocate (x(n)); call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Toeplitz matrix-vector product failed.")
         if (allocated(error)) return
      end block

      ! Matrix-matrix product.
      block
         real(dp), allocatable :: x(:, :), y(:, :), y_dense(:, :)
         allocate (x(n, n))
         call random_number(x)
         y = matmul(A, x); y_dense = matmul(dense(A), x)
         call check(error, all_close(y, y_dense), &
                    "Toeplitz matrix-matrix product failed.")
      end block
      return
   end subroutine test_matmul

   ! subroutine test_solve(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:)
   !
   !    ! Initialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !
   !    ! Solve with a single right-hand side vector.
   !    block
   !       real(dp), allocatable :: x(:), x_stdlib(:), b(:)
   !       allocate (b(n))
   !       ! Random rhs.
   !       call random_number(b)
   !       ! Solve with SpecialMatrices.
   !       x = solve(A, b)
   !       ! Solve with stdlib (dense).
   !       x_stdlib = solve(dense(A), b)
   !       ! Check error.
   !       call check(error, all_close(x, x_stdlib), &
   !                  "toeplitz solve with a single rhs failed.")
   !       if (allocated(error)) return
   !    end block
   !
   !    ! Solve with multiple right-hand side vectors.
   !    block
   !       real(dp), allocatable :: x(:, :), x_stdlib(:, :), b(:, :)
   !       allocate (b(n, n))
   !       ! Random rhs.
   !       call random_number(b)
   !       ! Solve with SpecialMatrices.
   !       x = solve(A, b)
   !       ! Solve with stdlib (dense).
   !       x_stdlib = solve(dense(A), b)
   !       ! Check error.
   !       call check(error, all_close(x, x_stdlib, abs_tol=1e-12_dp), &
   !                  "toeplitz solve with multiple rhs failed.")
   !    end block
   !
   !    return
   ! end subroutine test_solve

   ! subroutine test_svdvals(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:)
   !    real(dp), allocatable :: s(:), s_stdlib(:)
   !
   !    ! Initialize array.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !    ! Compute singular values.
   !    s = svdvals(A); s_stdlib = svdvals(dense(A))
   !    ! Check error.
   !    call check(error, all_close(s, s_stdlib), &
   !               "toeplitz svdvals failed.")
   !    return
   ! end subroutine test_svdvals
   !
   ! subroutine test_svd(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:), Amat(:, :)
   !    real(dp), allocatable :: u(:, :), s(:), vt(:, :)
   !
   !    ! Initialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !    ! Compute singular value decomposition.
   !    call svd(A, u, s, vt)
   !    ! Check error.
   !    allocate (Amat(n, n)); Amat = 0.0_dp
   !    Amat = matmul(u, matmul(diag(s), vt))
   !    call check(error, all_close(dense(A), Amat), &
   !               "toeplitz svd failed.")
   !    return
   ! end subroutine test_svd
   !
   ! subroutine test_eigvals(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:), abs_lambda(:)
   !    complex(dp), allocatable :: lambda(:), lambda_stdlib(:)
   !    complex(dp), allocatable :: lambda_check(:)
   !    integer(ilp) :: i, j
   !
   !    ! Initialize array.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !    ! Compute singular values.
   !    lambda = eigvals(A); lambda_stdlib = eigvals(dense(A))
   !    ! Re-order eigenvalues for comparison.
   !    allocate(lambda_check(n)) ; lambda_check = 0.0_dp
   !    do i = 1, n
   !       do j = 1, n
   !          if (is_close(lambda(i), lambda_stdlib(j))) then
   !             lambda_check(j) = lambda_stdlib(j)
   !             exit
   !          endif
   !       enddo
   !    enddo
   !    ! Check error.
   !    call check(error, all_close(lambda_check, lambda_stdlib), &
   !               "toeplitz eigvals failed.")
   !    return
   ! end subroutine test_eigvals
   !
   ! subroutine test_eig(error)
   !    type(error_type), allocatable, intent(out) :: error
   !    type(toeplitz) :: A
   !    real(dp), allocatable :: dv(:), Amat(:, :)
   !    complex(dp), allocatable :: lambda(:), left(:, :), right(:, :)
   !
   !    ! Initialize matrix.
   !    allocate (dv(n)); call random_number(dv); A = toeplitz(dv)
   !    ! Compute singular value decomposition.
   !    allocate (lambda(n), left(n, n), right(n, n))
   !    call eig(A, lambda, left=left, right=right)
   !    ! Check error.
   !    allocate (Amat(n, n)); Amat = 0.0_dp
   !    Amat = real(matmul(right, matmul(diag(lambda), hermitian(left))))
   !    call check(error, all_close(dense(A), Amat), &
   !               "toeplitz eig failed.")
   !    return
   ! end subroutine test_eig

end module
