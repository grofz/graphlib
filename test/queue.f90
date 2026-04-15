! Unit tests for a "queue" module.
! Can be used for a homework assignment in "Adv.Prog.Course"
!
! Files providing these modules must be compiled before compiling this file:
! - utest_mod: Compile file utest.f90.
! - queue_mod: Compile your own file queue.f90 that provides the implementation
!              of the queue. If the provided skeleton is used, it should compile
!              and run without errors (but the tests will fail)
!
! To compile and run the test, for example,
! $ gfortran utest.f90 queue.f90 main_test.f90
! % ./a.out
!
! Please read "README.txt".
! You do not have to read or modify this file or utest.f90 !
!

  program main
    use utest_mod
    use iso_fortran_env, only : dp => real64
    implicit none (type, external)
    type(utest_t) :: test
    integer :: maxsize, i
    integer, parameter :: MAXSIZE_SMALL = 1013, MAXSIZE_LARGE = 9800301
    integer, parameter :: ASK_REPEAT = 5
    real(dp) :: tstart, tend

    ! If copy constructor is not implemented, tests 21-23 must be skipped
    logical :: copyconstruct_implemented

    ! Not using module to make compilation easier
    interface
      subroutine run_test_queue(t, itest, maxsize)
        use conts_mod, only : queue_t
        import utest_t
        implicit none
        type(utest_t), intent(inout) :: t
        integer, intent(in) :: itest, maxsize
      end subroutine
    end interface

    ! Ask user
    print '("Ready to run queue tests?")'
    print '("Enter 1 if copy-constructor has been implemented and additional tests can run.")'
    print '("Otherwise, enter 0 to run the basic tests.")'
    QUESTION: block
      integer :: ios, choice
      do i = 1, ASK_REPEAT
        write(*,'("Enter 0 or 1: ")',advance='no')
        read(*,*,iostat=ios) choice
        if (ios==0) then
          select case(choice)
          case(1)
            copyconstruct_implemented = .true.
            exit QUESTION
          case(0)
            copyconstruct_implemented = .false.
            exit QUESTION
          case default
            continue
          end select
        end if
      end do
      stop 1 ! user was not able provide valid answer
    end block QUESTION
   
    ! Ask another question
    print '(/,"Ready for stress test with a large number of items?")'
    print '("Enter 1 to run a quick test. (",i0," items)")', MAXSIZE_SMALL
    print '("Enter 2 to run a stress test. (",i0," items)")', MAXSIZE_LARGE
    QUESTION2: block
      integer :: ios, choice
      do i = 1, ASK_REPEAT
        write(*,'("Enter 1 or 2: ")',advance='no')
        read(*,*,iostat=ios) choice
        if (ios==0) then
          select case(choice)
          case(1)
            maxsize = MAXSIZE_SMALL
            exit QUESTION2
          case(2)
            maxsize = MAXSIZE_LARGE
            exit QUESTION2
          case default
            continue
          end select
        end if
      end do
      stop 1 ! user was not able provide valid answer
    end block QUESTION2

    ! Now run the tests
    call cpu_time(tstart)
    print '("Running tests...",/)'

    test = utest_t()
    do i= 1, 23
      if (i<21 .or. copyconstruct_implemented) then
        call run_test_queue(test, i, MAXSIZE)
      else
        print '("Test ",i0," skipped because copy-constructor is not ready")', i
      end if
    enddo
    call cpu_time(tend)
    call test % summarize()
    print '("Time elapsed ",f9.3," s.")', tend-tstart
  end program main


! Note: The following procedures should be written in a module and a
! separate file.
! They are left as external procedures just to make compilation easier.

  subroutine run_test_queue(t, itest, maxsize)
    use utest_mod
    use conts_mod, only : queue_t, ERR_OK, ERR_EMPTY
    implicit none
    type(utest_t), intent(inout) :: t
    integer, intent(in) :: itest, maxsize
!
! Run unit tests for "queue" module
!
    type(queue_t) :: aq, bq, cq
    integer :: ierr, i
    logical :: isok
    integer, allocatable :: out(:)

    print '(a,i0)', 'Running test: ', itest

   !call aq%initialize(chunksize=1, capacity=1)
    call aq%initialize(chunksize=1, capacity=maxsize)
    if (itest==23 .or. itest==22) then
      call cq%initialize(chunksize=1, capacity=maxsize)
    else
      call cq%initialize(chunksize=1, capacity=1)
    end if

    select case(itest)
    case(1) ! empty queue is empty
      call t % assert_eq(aq % empty(), .true., ' 1. empty queue is empty')!

    case(2) ! empty queue has zero elements
      call t % assert_eq(aq % size(), 0, ' 2. empty queue has zero elements')

    case(3) ! remove from empty queue fails
      ierr = 0
      out = aq % dequeue(ierr)
      call t % assert_eq(ierr==ERR_EMPTY, .true., ' 3. remove from empty queue fails')

    case(4) ! peek at empty queue fails
      ierr = 0
      out = aq % peek(ierr)
      call t % assert_eq(ierr==ERR_EMPTY, .true., ' 4. peek at empty queue fails')

    case(6) ! adding one element, queue is no longer empty
      call aq % enqueue([100], ierr)
      call t % assert_eq(aq % empty(), .false., &
      & ' 6. adding one element, queue is no longer empty')

    case(7) ! adding one element, queue size is "1"
      call aq % enqueue([100], ierr)
      call t % assert_eq(aq % size(), 1, &
      & ' 7. adding one element, queue size is "1"')

    case(8) ! adding one element, peek works
      ierr = -1
      call aq % enqueue([100], ierr)
      out = aq%peek()
      isok = .false.
      if (allocated(out)) then
        if (size(out)==1) then
          isok = .true.
          call t % assert_eq(out(1), 100, &
          & ' 8. peek returns same number as added')
        end if
      end if
      call t % assert_eq(isok, .true., ' 8. peek returns correct size')
      call t % assert_eq(ierr, ERR_OK, ' 8. peek does give OK flag')
      call t % assert_eq(aq % size(), 1, ' 8. peek does not change size')

    case(9) ! adding one element, remove works
      call aq % enqueue([100], ierr)
      out = aq % dequeue(ierr)
      isok = .false.
      if (allocated(out)) then
        if (size(out)==1) then
          isok = .true.
          call t % assert_eq(out(1), 100, &
          & ' 9. remove returns expected item')
        end if
      end if
      call t % assert_eq(isok, .true., ' 9. remove returns correct size')
      call t % assert_eq(ierr, ERR_OK, ' 9. remove gives OK flag')
      call t % assert_eq(aq % empty(), .true., &
      & ' 9. after remove, queue is empty')
      call t % assert_eq(aq % size(), 0, ' 9. after remove, size is "0"')
      ierr = 0
      out = aq % dequeue(ierr)
      call t % assert_eq(ierr==ERR_EMPTY, .true., ' 9. repeated remove fails')

    case(10) ! adding maxsize elements works
      call add_items(aq, maxsize, isok)
      call t % assert_eq(isok, .true., '10. adding maximum elements works')
      call t % assert_eq(aq % size(), maxsize, &
      & '10. adding maximum elements, size is ok')

    case(12) ! adding maxsize elements, peek works
      call add_items(aq, maxsize, isok)
      out = aq % peek(ierr)
      call t % assert_eq(out(1), 1, &
      & '12. adding maximum elements, peek returns expected number')
      call t % assert_eq(ierr, ERR_OK, '12. peek gives ERR_OK flag')
      call t % assert_eq(aq % size(), maxsize, &
      & '12. peek does not change size')

    case(13) ! adding maxsize and removing all
      call add_items(aq, maxsize, isok)
      call remove_items(aq, 1, maxsize, isok)
      call t % assert_eq(isok, .true., '13. removing from full queue works')
      call t % assert_eq(aq % empty(), .true., '13. queue is empty again')

    case(14) ! adding full, removing half, adding third, removing all but one
      call add_items(aq, maxsize, isok)
      call remove_items(aq, 1, maxsize/2, isok)
      call t % assert_eq(isok, .true., '14. add full, remove half works')
      call add_items(aq, maxsize/3, isok)
      call remove_items(aq, maxsize/2+1, maxsize, isok)
      call t % assert_eq(isok, .true., '14. and add third, remove half works')
      call remove_items(aq, 1, maxsize/3-1, isok)
      call t % assert_eq(isok, .true., '14. and remove rest but one works')
      call t % assert_eq(aq % size(), 1, '14. one element remains')
      out = aq%peek(ierr)
      call t % assert_eq(out(1), maxsize/3, &
      & '14. remaining element correct')

    case(15) ! flow one element queue
      call aq % enqueue([1], ierr)
      call flow_queue(aq, 1, 1, 3*maxsize, isok)
      call t % assert_eq(isok, .true., '15. flow one element queue ok')

    case(16) ! flow halfsize queue
      call add_items(aq, maxsize/2, isok)
      call flow_queue(aq, 1, maxsize/2, 3*maxsize, isok)
      call t % assert_eq(isok, .true., '16. flow half queue ok')

    case(17) ! flow fullsize queue
      call add_items(aq, maxsize, isok)
      call flow_queue(aq, 1, maxsize, 3*maxsize, isok)
      call t % assert_eq(isok, .true., '17. flow full queue ok')

    case(21) ! copy constructor - copy is identical to original
      call add_items(aq, maxsize/2, isok)
      bq = aq
      call equal_queues(aq, bq, isok)
      call t % assert_eq(isok, .true., '21. copy is same as original')

    case(22) ! copy constructor - modifying original keeps copy intact
      call add_items(aq, maxsize/2, isok)
      call add_items(cq, maxsize/2, isok)
      bq = aq
      call remove_items(aq, 1, maxsize/3, isok)
      call equal_queues(bq, cq, isok)
      call t % assert_eq(isok, .true., &
      & '22. modifying original keeps copy intact')

    case(23) ! copy constructor - modifying copy keeps original intact
      call add_items(aq, maxsize/2, isok)
      call add_items(cq, maxsize/2, isok)
      bq = aq
      call remove_items(bq, 1, maxsize/3, isok)
      call equal_queues(aq, cq, isok)
      call t % assert_eq(isok, .true., &
      & '23. modifying copy keeps original intact')

    case(5,11,18,19,20) ! these tests are no longer used
      print '("Test ",i0," is no longer used")',itest

    case default
      print *, 'Warning: Test ',itest,' does not exist.'
    end select

  end subroutine run_test_queue



  subroutine add_items(q, n, isok)
    use conts_mod, only : queue_t
    implicit none
    type(queue_t), intent(inout) :: q
    integer, intent(in) :: n
    logical, intent(out) :: isok
!
! Add "n" items to the queue.
!
    integer :: i, ierr

    isok = .true.
    do i = 1, n
      call q % enqueue([i], ierr)
      if (ierr == 0) cycle
      isok = .false.
      exit
    enddo
  end subroutine add_items



  subroutine remove_items(q, i1, in, isok)
    use conts_mod, only : queue_t
    implicit none
    type(queue_t), intent(inout) :: q
    integer, intent(in) :: i1, in
    logical, intent(out) :: isok
!
! Test if expected items are being removed from the queue.
!
    integer :: i, ierr, nexpected
    integer :: irem(1)

    isok = .true.
    nexpected = q%size()
    do i = i1, in
      irem = q % dequeue(ierr)
      if (ierr == 0) nexpected = nexpected - 1
      if (ierr == 0 .and. irem(1) == i .and. q%size() == nexpected) cycle
      isok = .false.
      print '(a,i0,a,i0,a,i0,a,i0)', 'Remove fails: irem=',irem(1), &
      & ' expected=',i,' n=',q%size(),' ierr=',ierr
      exit
    enddo
  end subroutine remove_items



  subroutine equal_queues(q1, q2, isequal)
    use conts_mod, only : queue_t, ERR_OK
    implicit none
    type(queue_t), intent(inout) :: q1, q2
    logical, intent(out) :: isequal
!
! Test if two queues are identical. Note, the queues will be destroyed.
!
    integer :: i1(1), i2(1), ierr1, ierr2, safe_counter
    integer, parameter :: SAFE_EXIT = 10000000

    isequal = .false.
    if (q1 % size() > SAFE_EXIT) then
      print '(a,i0)', 'Equal_queues: this test fails because I can''t work with queues of size larger than ',SAFE_EXIT
      print '(a)', 'increase SAFE_EXIT or decrease queue size and try again'
      print *
      return
    endif
    safe_counter = 0 ! to avoid infinite loop

    if (q1 % size() /= q2 % size()) return
    do
      if (q1 % size() /= q2 % size()) exit
      if (q1 % empty() .and. q2 % empty()) then
        isequal = .true.
        exit
      endif
      if (q1 % empty() .or. q2 % empty()) exit
      i1 = q1 % dequeue(ierr1)
      i2 = q2 % dequeue(ierr2)
      if (any(i1 /= i2) .or. ierr1 /= ERR_OK .or. ierr2 /= ERR_OK) exit
      if (safe_counter > SAFE_EXIT) then
        print *, 'Equal_queues: infinite loop suspected - emeregency exit called'
        exit
      endif
      safe_counter = safe_counter + 1
    enddo
  end subroutine equal_queues



  subroutine flow_queue(q, i1, in, ncycles, isok)
    use conts_mod, only : queue_t, ERR_OK
    implicit none
    type(queue_t), intent(inout) :: q
    integer, intent(in) :: i1, in, ncycles
    logical, intent(out) :: isok
!
! Repeatedly remove and add elements to the queue
!
    integer :: ifirst, ilast, i, n, irem(1), ierr

    ifirst = i1
    ilast = in
    n = q % size()
    isok = .true.
    do i=1, ncycles
      irem = q % dequeue(ierr)
      if (ierr /= ERR_OK .or. irem(1) /= ifirst) then
        isok = .false.
        exit
      endif
      ifirst = ifirst + 1
      ilast = ilast + 1
      call q % enqueue([ilast], ierr)
      if (ierr /= 0) then
        isok = .false.
        exit
      endif
    enddo
  end subroutine flow_queue
