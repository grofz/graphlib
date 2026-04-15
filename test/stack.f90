! Main program is at the end
!
! Unit tests for a "stack_t" object.
!
  module stack_test_mod
    use utest_mod, only : utest_t
    use conts_mod, only : stack_t, ERR_OK, ERR_EMPTY
    implicit none
    private
    public run_test

  contains

    subroutine run_test(t, itest, maxsize)
      type(utest_t), intent(inout) :: t
      integer, intent(in) :: itest, maxsize
!
! Run unit tests for "stack_t" objects
!
      type(stack_t) :: aa, bb
      integer :: ierr, i, j
      integer, allocatable :: arr(:), out(:)
      logical :: isok

      print '(a,i0)', 'Running test: ', itest

      ! Initialization of stack?
      call aa%initialize(chunksize=1, capacity=1)
      ierr = 0
      i = 0

      if (itest > 5) then
        ! Add many items to stack
        allocate(arr(maxsize))
        do i=1,maxsize
          arr(i) = i
        end do
        do i=1,maxsize
          call aa % push(arr(i:i))
        end do
        isok = .true.
      end if


      ! Select test
      select case(itest)
      case(1) ! empty stack is empty
        call t % assert_eq(aa%empty(), .true., ' 1. new stack is empty')

      case(2) ! putting element to stack: stack no longer empty
        call aa % push([10])
        call t % assert_eq(aa%empty(), .false., ' 2. stack no longer empty')

      case(3) ! poping empty stack fails
        out = aa % pop(ierr)
        call t % assert_eq(ierr==ERR_EMPTY, .true., ' 3. poping empty stack fails')

      case(4) ! peeking empty stack fails
        out = aa % peek(ierr)
        call t % assert_eq(ierr==ERR_EMPTY, .true., ' 4. peeking empty stack fails')

      case(5) ! puting single item to stack and removing it
        call aa % push([10])
        i = -1
        ierr = -1
        allocate(out(0))
        out = aa % peek(ierr)
        call t % assert_eq(ierr==ERR_OK, .true., ' 5a. peeking stack no error')
        if (size(out)==1) then
          call t % assert_eq(out(1)==10, .true., ' 5b. peeking stack correct value')
        else
          call t % assert_eq(.false., .true., ' 5b. wrong dimension')
        end if
        deallocate(out)

        i = -1
        ierr = -1
        out = aa % pop(ierr)
        call t % assert_eq(ierr==ERR_OK, .true., ' 5c. poping stack no error')
        if (size(out)==1) then
          call t % assert_eq(out(1)==10, .true., ' 5d. poping stack correct value')
        else
          call t % assert_eq(.false., .true., ' 5c. wrong dimension')
        end if
        deallocate(out)
        call t % assert_eq(aa%empty(), .true., ' 5e. stack is empty again')

      case(6) ! puting many items to stack and removing in correct order
        do i=maxsize,maxsize/2,-1
          out = aa % peek(ierr)
          if (out(1) /= arr(i)) isok = .false.
          if (allocated(out)) deallocate(out)
          out = aa % pop(ierr)
          if (allocated(out)) then
            if (size(out)==1) then
              if (out(1)==arr(i)) cycle
            end if
          end if
          isok = .false.
        end do
        call t % assert_eq(isok, .true., ' 6. put/pop many elements ok')

      case(7) ! putting items / making copy / removing items / copy still ok
        bb = aa
        do i=maxsize,maxsize/2,-1
          out = aa % pop(ierr)
        end do
        do i=maxsize,1,-1
          out = bb%peek(ierr)
          if (out(1) /= arr(i)) isok = .false.
          out = bb % pop(ierr)
          if (out(1) /= arr(i)) isok = .false.
        end do
        call t % assert_eq(isok, .true., ' 7. copied stack ok')
        call t % assert_eq(bb%empty(), .true., ' 7b. stack empty')

      end select
    end subroutine run_test

  end module stack_test_mod



  program main
    call sub()
  contains
    subroutine sub()
      use utest_mod, only : utest_t
      use stack_test_mod, only : run_test
      implicit none
      type(utest_t) :: test
      integer :: i
      logical :: copy_construct_implemented
      !integer, parameter :: MAXSIZE = 10
      integer, parameter :: MAXSIZE = 1000000

      test = utest_t()

      print '("Ready to run stack tests?")'
      print '("Enter 1 if copy-constructor has been implemented and additional test can run.")'
      print '("Otherwise, enter 0 to safely run basic tests only.")'
      QUESTION: block
        integer :: ios, choice
        do i = 1, 10
          write(*,'("Enter 0 or 1: ")',advance='no')
          read(*,*,iostat=ios) choice
          if (ios==0) then
            select case(choice)
            case(1)
              copy_construct_implemented = .true.
              exit QUESTION
            case(0)
              copy_construct_implemented = .false.
              exit QUESTION
            case default
              continue
            end select
          end if
        end do
        stop 1 ! user was not able provide valid answer
      end block QUESTION

      ! Do not run the 7th test until assignment is defined
      do i= 1, 7
        if (i==7 .and. .not. copy_construct_implemented) then
          print '("Test 7 skiped")'
          cycle
        end if
        call run_test(test, i, MAXSIZE)
      enddo
      call test % summarize()
    end subroutine sub
  end program main
