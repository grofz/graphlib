  module conts_test_mod
    use utest_mod, only : utest_t
    use conts_mod, only : stack_t, queue_t
    implicit none (type, external)
    private
    public run_stack_test, run_queue_test
  contains

    subroutine run_stack_test(t, itest, maxsize)
      type(utest_t), intent(inout) :: t
      integer, intent(in) :: itest, maxsize
!
! Run unit tests for "stack_t" objects
!
      type(stack_t) :: aa, bb
      integer :: i
      integer, allocatable :: arr(:), out(:)
      logical :: isok

      print '(a,i0)', 'Running test: ', itest

      ! Initialization of stack?
      call aa%initialize(chunksize=1, capacity=1)
      i = 0

      if (itest > 3) then
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
        call t % assert(aa%empty(), .true., ' 1. new stack is empty')

      case(2) ! putting element to stack: stack no longer empty
        call aa % push([10])
        call t % assert(aa%empty(), .false., ' 2. stack no longer empty')

      case(3) ! puting single item to stack and removing it
        call aa % push([10])
        i = -1
        allocate(out(0))
        out = aa % peek()
        if (size(out)==1) then
          call t % assert(out(1)==10, .true., ' 3a. peeking stack correct value')
        else
          call t % assert(.false., .true., ' 3a. wrong dimension')
        end if
        deallocate(out)

        i = -1
        out = aa % pop()
        if (size(out)==1) then
          call t % assert(out(1)==10, .true., ' 3b. poping stack correct value')
        else
          call t % assert(.false., .true., ' 3b. wrong dimension')
        end if
        deallocate(out)
        call t % assert(aa%empty(), .true., ' 3c. stack is empty again')

      case(4) ! puting many items to stack and removing in correct order
        do i=maxsize,maxsize/2,-1
          out = aa % peek()
          if (out(1) /= arr(i)) isok = .false.
          if (allocated(out)) deallocate(out)
          out = aa % pop()
          if (allocated(out)) then
            if (size(out)==1) then
              if (out(1)==arr(i)) cycle
            end if
          end if
          isok = .false.
        end do
        call t % assert(isok, .true., ' 4. put/pop many elements ok')

      case(5) ! putting items / making copy / removing items / copy still ok
        bb = aa
        do i=maxsize,maxsize/2,-1
          out = aa % pop()
        end do
        do i=maxsize,1,-1
          out = bb%peek()
          if (out(1) /= arr(i)) isok = .false.
          out = bb % pop()
          if (out(1) /= arr(i)) isok = .false.
        end do
        call t % assert(isok, .true., ' 5. copied stack ok')
        call t % assert(bb%empty(), .true., ' 5b. stack empty')

      end select
    end subroutine run_stack_test


    subroutine run_queue_test(t, itest, maxsize)
      type(utest_t), intent(inout) :: t
      integer, intent(in) :: itest, maxsize
!
! Run unit tests for "queue" module
!
      type(queue_t) :: aq, bq, cq
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
        call t % assert(aq % empty(), .true., ' 1. empty queue is empty')!

      case(2) ! empty queue has zero elements
        call t % assert(aq % size(), 0, ' 2. empty queue has zero elements')

      case(3) ! remove from empty queue fails
        continue ! not relevant anymore

      case(4) ! peek at empty queue fails
        continue ! not relevant anymore

      case(6) ! adding one element, queue is no longer empty
        call aq % enqueue([100])
        call t % assert(aq % empty(), .false., &
        & ' 6. adding one element, queue is no longer empty')

      case(7) ! adding one element, queue size is "1"
        call aq % enqueue([100])
        call t % assert(aq % size(), 1, &
        & ' 7. adding one element, queue size is "1"')

      case(8) ! adding one element, peek works
        call aq % enqueue([100])
        out = aq%peek()
        isok = .false.
        if (allocated(out)) then
          if (size(out)==1) then
            isok = .true.
            call t % assert(out(1), 100, &
            & ' 8. peek returns same number as added')
          end if
        end if
        call t % assert(isok, .true., ' 8. peek returns correct size')
        call t % assert(aq % size(), 1, ' 8. peek does not change size')

      case(9) ! adding one element, remove works
        call aq % enqueue([100])
        out = aq % dequeue()
        isok = .false.
        if (allocated(out)) then
          if (size(out)==1) then
            isok = .true.
            call t % assert(out(1), 100, &
            & ' 9. remove returns expected item')
          end if
        end if
        call t % assert(isok, .true., ' 9. remove returns correct size')
        call t % assert(aq % empty(), .true., &
        & ' 9. after remove, queue is empty')
        call t % assert(aq % size(), 0, ' 9. after remove, size is "0"')

      case(10) ! adding maxsize elements works
        call add_items(aq, maxsize, isok)
        call t % assert(isok, .true., '10. adding maximum elements works')
        call t % assert(aq % size(), maxsize, &
        & '10. adding maximum elements, size is ok')

      case(12) ! adding maxsize elements, peek works
        call add_items(aq, maxsize, isok)
        out = aq % peek()
        call t % assert(out(1), 1, &
        & '12. adding maximum elements, peek returns expected number')
        call t % assert(aq % size(), maxsize, &
        & '12. peek does not change size')

      case(13) ! adding maxsize and removing all
        call add_items(aq, maxsize, isok)
        call remove_items(aq, 1, maxsize, isok)
        call t % assert(isok, .true., '13. removing from full queue works')
        call t % assert(aq % empty(), .true., '13. queue is empty again')

      case(14) ! adding full, removing half, adding third, removing all but one
        call add_items(aq, maxsize, isok)
        call remove_items(aq, 1, maxsize/2, isok)
        call t % assert(isok, .true., '14. add full, remove half works')
        call add_items(aq, maxsize/3, isok)
        call remove_items(aq, maxsize/2+1, maxsize, isok)
        call t % assert(isok, .true., '14. and add third, remove half works')
        call remove_items(aq, 1, maxsize/3-1, isok)
        call t % assert(isok, .true., '14. and remove rest but one works')
        call t % assert(aq % size(), 1, '14. one element remains')
        out = aq%peek()
        call t % assert(out(1), maxsize/3, &
        & '14. remaining element correct')

      case(15) ! flow one element queue
        call aq % enqueue([1])
        call flow_queue(aq, 1, 1, 3*maxsize, isok)
        call t % assert(isok, .true., '15. flow one element queue ok')

      case(16) ! flow halfsize queue
        call add_items(aq, maxsize/2, isok)
        call flow_queue(aq, 1, maxsize/2, 3*maxsize, isok)
        call t % assert(isok, .true., '16. flow half queue ok')

      case(17) ! flow fullsize queue
        call add_items(aq, maxsize, isok)
        call flow_queue(aq, 1, maxsize, 3*maxsize, isok)
        call t % assert(isok, .true., '17. flow full queue ok')

      case(21) ! copy constructor - copy is identical to original
        call add_items(aq, maxsize/2, isok)
        bq = aq
        call equal_queues(aq, bq, isok)
        call t % assert(isok, .true., '21. copy is same as original')

      case(22) ! copy constructor - modifying original keeps copy intact
        call add_items(aq, maxsize/2, isok)
        call add_items(cq, maxsize/2, isok)
        bq = aq
        call remove_items(aq, 1, maxsize/3, isok)
        call equal_queues(bq, cq, isok)
        call t % assert(isok, .true., &
        & '22. modifying original keeps copy intact')

      case(23) ! copy constructor - modifying copy keeps original intact
        call add_items(aq, maxsize/2, isok)
        call add_items(cq, maxsize/2, isok)
        bq = aq
        call remove_items(bq, 1, maxsize/3, isok)
        call equal_queues(aq, cq, isok)
        call t % assert(isok, .true., &
        & '23. modifying copy keeps original intact')

      case(5,11,18,19,20) ! these tests are no longer used
        print '("Test ",i0," is no longer used")',itest

      case default
        print *, 'Warning: Test ',itest,' does not exist.'
      end select

    end subroutine run_queue_test


    subroutine add_items(q, n, isok)
      type(queue_t), intent(inout) :: q
      integer, intent(in) :: n
      logical, intent(out) :: isok
!
! Add "n" items to the queue.
!
      integer :: i

      isok = .true.
      do i = 1, n
        call q % enqueue([i])
      enddo
    end subroutine add_items


    subroutine remove_items(q, i1, in, isok)
      type(queue_t), intent(inout) :: q
      integer, intent(in) :: i1, in
      logical, intent(out) :: isok
!
! Test if expected items are being removed from the queue.
!
      integer :: i, nexpected
      integer :: irem(1)

      isok = .true.
      nexpected = q%size()
      do i = i1, in
        irem = q % dequeue()
        nexpected = nexpected - 1
        if (irem(1) == i .and. q%size() == nexpected) cycle
        isok = .false.
        print '(a,i0,a,i0,a,i0)', 'Remove fails: irem=',irem(1), &
        & ' expected=',i,' n=',q%size()
        exit
      enddo
    end subroutine remove_items


    subroutine equal_queues(q1, q2, isequal)
      type(queue_t), intent(inout) :: q1, q2
      logical, intent(out) :: isequal
!
! Test if two queues are identical. Note, the queues will be destroyed.
!
      integer :: i1(1), i2(1), safe_counter
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
        i1 = q1 % dequeue()
        i2 = q2 % dequeue()
        if (any(i1 /= i2)) exit
        if (safe_counter > SAFE_EXIT) then
          print *, 'Equal_queues: infinite loop suspected - emeregency exit called'
          exit
        endif
        safe_counter = safe_counter + 1
      enddo
    end subroutine equal_queues


    subroutine flow_queue(q, i1, in, ncycles, isok)
      type(queue_t), intent(inout) :: q
      integer, intent(in) :: i1, in, ncycles
      logical, intent(out) :: isok
!
! Repeatedly remove and add elements to the queue
!
      integer :: ifirst, ilast, i, n, irem(1)

      ifirst = i1
      ilast = in
      n = q % size()
      isok = .true.
      do i=1, ncycles
        irem = q % dequeue()
        if (irem(1) /= ifirst) then
          isok = .false.
          exit
        endif
        ifirst = ifirst + 1
        ilast = ilast + 1
        call q % enqueue([ilast])
      enddo
    end subroutine flow_queue

  end module conts_test_mod
