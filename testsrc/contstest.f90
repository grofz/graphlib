!
! Unit tests for stack and queue
!
  program main
    call sub()
  contains
    subroutine sub()
      use iso_fortran_env, only : dp=>real64
      use utest_mod, only : utest_t
      use conts_test_mod, only : run_stack_test, run_queue_test
      implicit none (type, external)
      type(utest_t) :: test_stack, test_queue
      integer :: i, maxsize
      integer, parameter :: MAXSIZE_SMALL = 1013, MAXSIZE_LARGE = 9800301
      real(dp) :: tstart, tend

      maxsize = MAXSIZE_LARGE
      test_stack = utest_t()
      test_queue = utest_t()
      call cpu_time(tstart)

      print '(/,"Runing stack tests...")'
      do i= 1, 5
        call run_stack_test(test_stack, i, MAXSIZE)
      enddo

      print '(/,"Running queue tests...")'
      do i= 1, 23
        call run_queue_test(test_queue, i, MAXSIZE)
      enddo

      call cpu_time(tend)
      print '(/,"Time elapsed ",f9.3," s.")', tend-tstart

      print '(/,"Stack results")'
      call test_stack % summarize()
      print '(/,"Queue results")'
      call test_queue % summarize()
    end subroutine sub

  end program main
