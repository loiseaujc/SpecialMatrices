program Tester
   ! For best practice.
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   ! Unit-test utility.
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   ! List of tests.
   use test_diagonal, only: collect_diagonal_testsuite
   use test_bidiagonal, only: collect_bidiagonal_testsuite
   use test_tridiagonal, only: collect_tridiagonal_testsuite
   use test_symtridiagonal, only: collect_symtridiagonal_testsuite
   use test_strang, only: collect_strang_testsuite
   use test_poisson2D, only: collect_poisson2D_testsuite
   use test_circulant, only: collect_circulant_testsuite
   use test_toeplitz, only: collect_toeplitz_testsuite
   use test_hankel, only: collect_hankel_testsuite

   implicit none(type, external)

   ! Unit-test related.
   integer :: status, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   status = 0

   testsuites = [ &
                new_testsuite("Diagonal Matrices", collect_diagonal_testsuite), &
                new_testsuite("Bidiagonal Matrices", collect_bidiagonal_testsuite), &
                new_testsuite("Tridiagonal Matrices", collect_tridiagonal_testsuite), &
                new_testsuite("SymTridiagonal Matrices", collect_symtridiagonal_testsuite), &
                new_testsuite("Strang Matrices", collect_strang_testsuite), &
                new_testsuite("Poisson2D Matrices", collect_poisson2D_testsuite), &
                new_testsuite("Circulant Matrices", collect_circulant_testsuite), &
                new_testsuite("Toeplitz Matrices", collect_toeplitz_testsuite), &
                new_testsuite("Hankel Matrices", collect_hankel_testsuite) &
                ]

   do is = 1, size(testsuites)
      write (output_unit, *) "-----"
      write (output_unit, fmt) "Testing :", testsuites(is)%name
      write (output_unit, *) "-----", new_line("a")
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (output_unit, *) new_line("a")
   end do

   if (status > 0) then
      write (error_unit, "(i0, 1x, a)") status, "test(s) failed!"
   else if (status == 0) then
      write (output_unit, *) "All tests succesfully passed!", new_line("a")
   end if

end program Tester
