module conts_mod
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none (type, external)
  private

  type, public :: stack_t
    private
    integer :: chunksize = 1
    integer :: capacity = 0
    integer :: n = 0
    integer, allocatable :: arr(:)
  contains
    procedure :: initialize => stack_initialize
    procedure :: push => stack_push, pop => stack_pop
    procedure :: size => stack_size, empty => stack_empty
  end type
  ! TODO: clear


  type, public :: queue_t
    private
    integer :: chunksize = 1
    integer :: capacity = 0
    integer :: n = 0
    integer :: rear
    integer, allocatable :: arr(:)
  contains
    procedure :: initialize => queue_initialize
    procedure :: enqueue => queue_enqueue, dequeue => queue_dequeue
    procedure :: size => queue_size, empty => queue_empty
  end type
  ! TODO: clear


  integer, parameter :: HMAP_NULL = -1
  integer, parameter :: PQUEUE_MIN = 1, & ! lower P is higher priority
                        PQUEUE_MAX = 2    ! higher P is higher priority
  integer, parameter :: INTEGER_MOLD(0) = [integer ::]

  type, public :: handle_t
    private
    integer :: index_to_hmap = HMAP_NULL
    integer :: version = 1
  contains
    procedure, private :: handle_eq
    generic :: operator(==) => handle_eq
  end type

  type, public :: pqueue_t
    integer, allocatable :: values(:,:)       ! heap of values
    integer, allocatable :: priorities(:)     ! heap of priorities
    type(handle_t), allocatable :: handles(:) ! heap of handles
    integer, allocatable :: hmap(:)           ! handle to map to item's position in heaps
    type(queue_t) :: free_handles
    integer :: n = 0
    integer :: chunksize = 1
    integer :: ordering = PQUEUE_MIN
  contains
    procedure :: initialize => pqueue_initialize
    procedure :: insert => pqueue_insert, pop => pqueue_pop
    procedure :: update_priority => pqueue_update_priority
    procedure :: size => pqueue_size, empty => pqueue_empty
    procedure :: contains => pqueue_contains
  end type
  ! TODO: remove(handle)
  ! TODO: clear
  ! TODO: get_priority(handle)


  integer :: efid = error_unit
  integer, parameter :: DEFAULT_CAPACITY = 10
  integer, parameter, public :: &
    ERR_OK = 0, &
    ERR_EMPTY = 1, &
    ERR_INVALID_HANDLE = 2, &
    ERR_INVALID_ARG_SIZE = 3

contains

! ----------------------------
! STACK / QUEUE implementation
! ----------------------------

  subroutine stack_initialize(this, chunksize, capacity)
    class(stack_t), intent(inout) :: this
    integer, intent(in), optional :: chunksize, capacity

    integer :: new_capacity

    if (allocated(this%arr)) deallocate(this%arr)
    if (present(chunksize)) then
      this%chunksize = chunksize
    else
      this%chunksize = 1
    end if
    if (present(capacity)) then
      new_capacity = max(1, capacity)
    else
      new_capacity = DEFAULT_CAPACITY
    end if
    this%n = 0
    allocate(this%arr(0))
    call stack_increase_capacity(this, new_capacity)
  end subroutine


  subroutine queue_initialize(this, chunksize, capacity)
    class(queue_t), intent(inout) :: this
    integer, intent(in), optional :: chunksize, capacity

    integer :: new_capacity

    if (allocated(this%arr)) deallocate(this%arr)
    if (present(chunksize)) then
      this%chunksize = chunksize
    else
      this%chunksize = 1
    end if
    if (present(capacity)) then
      new_capacity = max(1, capacity)
    else
      new_capacity = DEFAULT_CAPACITY
    end if
    this%n = 0
    this%rear = 1
    this%capacity = new_capacity
    allocate(this%arr(0))
    call queue_increase_capacity(this, new_capacity)
  end subroutine


  subroutine stack_increase_capacity(this, new_capacity)
    class(stack_t), intent(inout) :: this
    integer, intent(in) :: new_capacity

    integer, allocatable :: tmp(:)

    if (new_capacity <= this%n) &
      & error stop 'stack_increase_capacity - internal error'
    allocate(tmp(new_capacity*this%chunksize))
    tmp(1:this%n*this%chunksize) = this%arr(1:this%n*this%chunksize)
    call move_alloc(tmp, this%arr)
    this%capacity = new_capacity
print *, 'stack: new capacity ', new_capacity
  end subroutine stack_increase_capacity


  subroutine queue_increase_capacity(this, new_capacity)
    class(queue_t), intent(inout) :: this
    integer, intent(in) :: new_capacity

    integer, allocatable :: tmp(:)
    integer :: front, c1, c2

    if (new_capacity <= this%n) &
      & error stop 'queue_increase_capacity - internal error'
    allocate(tmp(new_capacity*this%chunksize))
    front = modulo(this%rear-this%n*this%chunksize-1, this%capacity*this%chunksize) + 1
    if (front >= this%rear .and. this%n/=0) then
      c1 = this%capacity*this%chunksize-front+1
      c2 = this%n*this%chunksize-c1
      tmp(1:c1) = this%arr(front:this%capacity*this%chunksize)
      tmp(c1+1:c1+c2) = this%arr(1:c2)
    else
      tmp(1:this%n*this%chunksize) = this%arr(front:this%rear-1)
    end if
    call move_alloc(tmp, this%arr)
    this%capacity = new_capacity
    this%rear = this%n*this%chunksize+1
print *, 'queue: new capacity ', new_capacity
  end subroutine queue_increase_capacity

!
! Notes: Circular array for the queue
! F = (R-N*CS-1 % CAP) + 1
!
! |--- --- --- --- --- ---| empty
! |R=F
!
! |### --- --- --- --- ---| enqueue 1
! |F   R                    F = 4-3-1 % 18 + 1
!
! |### ### --- --- --- ---| enqueue 2
! |F       R                F = 7-6-1 % 18 + 1
!
! |--- --- --- ### ### ###|
! |R           F            F = 1-9-1 % 18 + 1 = 10
!
! |### ### ### --- --- ###|
!              R       F    F = 10-12-1 % 18 + 1 = 16
!
! |### ### ### --- --- ---|
!  F           R            F = 10-9-1 % 18 + 1 = 1


  subroutine stack_push(this, newitem, ierr)
    class(stack_t), intent(inout) :: this
    integer, intent(in) :: newitem(:)
    integer, intent(out), optional :: ierr

    if (size(newitem) /= this%chunksize) then
      call handle_error(ERR_INVALID_ARG_SIZE, 'stack_push - unexpected size of newitem', ierr)
      return
    else if (this%capacity == 0) then
      error stop 'stack_push - stack is not initialized'
    else
      if (present(ierr)) ierr = ERR_OK
    end if
    if (this%capacity == this%n) &
      & call stack_increase_capacity(this, 2*this%capacity)

    ! push new item
    block
      integer :: is, ie
      this%n = this%n + 1
      is = (this%n-1)*this%chunksize+1
      ie = is + this%chunksize-1
      this%arr(is:ie) = newitem
    end block
  end subroutine stack_push


  subroutine queue_enqueue(this, newitem, ierr)
    class(queue_t), intent(inout) :: this
    integer, intent(in) :: newitem(:)
    integer, intent(out), optional :: ierr

    if (size(newitem) /= this%chunksize) then
      call handle_error(ERR_INVALID_ARG_SIZE, 'queue_enqueue - unexpected size of newitem', ierr)
      return
    else if (this%capacity == 0) then
      error stop 'queue_enqueue - queue is not itialized'
    else
      if (present(ierr)) ierr = ERR_OK
    end if
    if (this%capacity == this%n) &
      & call queue_increase_capacity(this, 2*this%capacity)

    ! enqueue new item
    block
      integer :: is, ie
      is = this%rear
      ie = is + this%chunksize-1
      if (ie<=is) error stop 'queue_enqueue - logical error'
      this%arr(is:ie) = newitem
      this%n = this%n + 1
      this%rear = modulo(this%rear+this%chunksize-1, this%capacity*this%chunksize)+1
    end block
  end subroutine queue_enqueue


  function stack_pop(this, ierr) result(pop_item)
    class(stack_t), intent(inout) :: this
    integer, intent(out), optional :: ierr
    integer :: pop_item(this%chunksize)

    if (this%n < 1) then
      call handle_error(ERR_EMPTY, 'stack_pop - empty stack', ierr)
      return
    else
      if (present(ierr)) ierr = ERR_OK
    end if

    block
      integer :: is, ie
      is = (this%n-1)*this%chunksize+1
      ie = is + this%chunksize-1
      pop_item = this%arr(is:ie)
      this%n = this%n - 1
    end block
  end function stack_pop


  function queue_dequeue(this, ierr) result(pop_item)
    class(queue_t), intent(inout) :: this
    integer, intent(out), optional :: ierr
    integer :: pop_item(this%chunksize)

    if (this%n < 1) then
      call handle_error(ERR_EMPTY, 'queue_dequeue - empty queue', ierr)
      return
    else
      if (present(ierr)) ierr = ERR_OK
    end if

    block
      integer :: is, ie
      is = modulo(this%rear-this%n*this%chunksize-1, this%capacity*this%chunksize) + 1
      ie = is + this%chunksize-1
      pop_item = this%arr(is:ie)
      this%n = this%n - 1
    end block
  end function queue_dequeue


  pure function stack_size(this) result(n)
    class(stack_t), intent(in) :: this
    integer :: n
    n = this%n
  end function stack_size


  pure function queue_size(this) result(n)
    class(queue_t), intent(in) :: this
    integer :: n
    n = this%n
  end function queue_size


  pure function stack_empty(this) result(is_empty)
    class(stack_t), intent(in) :: this
    logical :: is_empty
    is_empty = this%n == 0
  end function stack_empty


  pure function queue_empty(this) result(is_empty)
    class(queue_t), intent(in) :: this
    logical :: is_empty
    is_empty = this%n == 0
  end function queue_empty


  ! ------------------------------
  ! PRIORITY QUEUE implementation
  ! ------------------------------

  pure function handle_eq(a, b) result(eq)
    class(handle_t), intent(in) :: a, b
    logical :: eq
    eq = a%version==b%version .and. a%index_to_hmap==b%index_to_hmap
  end function handle_eq


  pure function is_higher_priority(pa, pb, this) result(is)
    integer, intent(in) :: pa, pb
    class(pqueue_t), intent(in) :: this
    logical :: is
    select case(this%ordering)
    case(PQUEUE_MIN) ! assuming lower "P" means higher priority...
      is = pa < pb
    case(PQUEUE_MAX) ! assuming higher "P" means higher priority...
      is = pa > pb
    case default
      error stop 'invalid ORDERING flag in pqueue'
    end select
  end function


  pure function is_lower_priority(pa, pb, this) result(is)
    integer, intent(in) :: pa, pb
    class(pqueue_t), intent(in) :: this
    logical :: is
    select case(this%ordering)
    case(PQUEUE_MIN) ! assuming lower "P" means higher priority...
      is = pa > pb
    case(PQUEUE_MAX) ! assuming higher "P" means higher priority...
      is = pa < pb
    case default
      error stop 'invalid ORDERING flag in pqueue'
    end select
  end function


  subroutine pqueue_initialize(this, chunksize, capacity, ordering)
    class(pqueue_t), intent(inout) :: this
    integer, intent(in), optional :: chunksize, capacity, ordering

    integer :: new_capacity
    type(handle_t) :: handle

    if (present(chunksize)) then
      this%chunksize = chunksize
    else
      this%chunksize = 1
    end if

    this%ordering = PQUEUE_MIN
    if (present(ordering)) then
      if (ordering==PQUEUE_MAX .or. ordering==PQUEUE_MIN) then
        this%ordering = ordering
      else
        error stop 'pqueue_initialize - invalid ordering tag'
      end if
    end if

    if (present(capacity)) then
      new_capacity = capacity
    else
      new_capacity = DEFAULT_CAPACITY
    end if

    if (allocated(this%values)) deallocate(this%values)
    allocate(this%values(chunksize,0))
    if (allocated(this%priorities)) deallocate(this%priorities)
    allocate(this%priorities(0))
    if (allocated(this%handles)) deallocate(this%handles)
    allocate(this%handles(0))
    if (allocated(this%hmap)) deallocate(this%hmap)
    allocate(this%hmap(0))

    call this%free_handles%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
    this%n = 0
    call pqueue_increase_capacity(this, new_capacity)
  end subroutine pqueue_initialize


  subroutine pqueue_increase_capacity(this, new_capacity)
    class(pqueue_t), intent(inout) :: this
    integer, intent(in), optional :: new_capacity

    integer :: old_capacity, new_capacity0
    integer, allocatable :: tmp2(:,:), tmp1(:)
    type(handle_t), allocatable :: tmp_handles(:)

    old_capacity = size(this%values, dim=2)
    if (present(new_capacity)) then
      new_capacity0 = new_capacity
    else
      new_capacity0 = 2*old_capacity
    end if

    allocate(tmp2(size(this%values,1), new_capacity0))
    tmp2(:,1:old_capacity) = this%values
    call move_alloc(tmp2, this%values)

    allocate(tmp1(new_capacity0))
    tmp1(1:old_capacity) = this%priorities
    call move_alloc(tmp1, this%priorities)

    allocate(tmp_handles(new_capacity0))
    tmp_handles(1:old_capacity) = this%handles
    call move_alloc(tmp_handles, this%handles)

    allocate(tmp1(new_capacity0), source=HMAP_NULL)
    tmp1(1:old_capacity) = this%hmap
    call move_alloc(tmp1, this%hmap)

    block ! create fresh handles
      integer :: i
      type(handle_t) :: new_handle
      do i=old_capacity+1, new_capacity0
        new_handle%index_to_hmap = i
        new_handle%version = 1
        call this%free_handles%enqueue(transfer(new_handle, INTEGER_MOLD))
      end do
    end block
  end subroutine pqueue_increase_capacity


  subroutine borrow_handle(this, handle)
    class(pqueue_t), intent(inout) :: this
    type(handle_t), intent(out) :: handle

    if (this%free_handles%size()==0) call pqueue_increase_capacity(this)
    if (this%free_handles%size()==0) error stop 'borrow_handle - no more handles available'
    handle = transfer(this%free_handles%dequeue(), handle)
  end subroutine borrow_handle


  subroutine return_handle(this, handle)
    class(pqueue_t), intent(inout) :: this
    type(handle_t), intent(in) :: handle

    type(handle_t) :: reused_handle

    reused_handle = handle
    reused_handle%version = reused_handle%version + 1
    call this%free_handles%enqueue(transfer(reused_handle,INTEGER_MOLD))
  end subroutine return_handle


  pure function get_idheap(this, handle) result(id)
    class(pqueue_t), intent(in) :: this
    type(handle_t), intent(in) :: handle
    integer :: id
!
! Return position of an item in queue using handle. If handle refers to an
! item that is no longer in queue, HMAP_NULL value is returned
!
    id = HMAP_NULL
    if (handle%index_to_hmap > 0 .and. handle%index_to_hmap <= size(this%hmap)) then
      id = this%hmap(handle%index_to_hmap)
      if (id/=HMAP_NULL) then
        ! verify version component matches the stored one
        if (this%handles(id)%version/=handle%version) id = HMAP_NULL
      end if
    end if
  end function get_idheap


  pure function pqueue_contains(this, handle) result(in_queue)
    class(pqueue_t), intent(in) :: this
    type(handle_t), intent(in) :: handle
    logical :: in_queue
!
! Use this function to test if "handle" references to an item still in
! the queue.
!
    in_queue = get_idheap(this, handle) /= HMAP_NULL
  end function pqueue_contains


  pure subroutine swap(this, i, j)
    type(pqueue_t), intent(inout) :: this
    integer, intent(in) :: i, j
!
! Swap heap elements ("values", "handles" and "priorities" arrays)
! at "i" and "j" positions and reflect this change in the hmap table.
!
    integer :: tmp(size(this%values,dim=1)), idi, idj
    type(handle_t) :: tmp_handle

    ! find positions of items in the hash-map
    idi = this%handles(i)%index_to_hmap
    idj = this%handles(j)%index_to_hmap

    ! swap values, handles and priorities
    tmp = this%values(:,i)
    this%values(:,i) = this%values(:,j)
    this%values(:,j) = tmp
    tmp_handle = this%handles(i)
    this%handles(i) = this%handles(j)
    this%handles(j) = tmp_handle
    tmp(1) = this%priorities(i)
    this%priorities(i) = this%priorities(j)
    this%priorities(j) = tmp(1)

    ! update hash-map
    this%hmap(idi) = j
    this%hmap(idj) = i
  end subroutine swap


  pure subroutine bubble_up(this, ind)
    type(pqueue_t), intent(inout) :: this
    integer, intent(in) :: ind
!
! Helper procedure to move items with higher priority up the heap.
!
    integer :: pari, curi

    pari = ind
    do while (pari > 1)
      curi = pari     ! current node
      pari = pari / 2 ! parent node
      if (is_lower_priority(this%priorities(pari), this%priorities(curi), this)) then
        ! swap nodes if parent node has lower priority than the current node
        call swap(this, curi, pari)
      else
        exit
      end if
    end do
  end subroutine bubble_up


  pure subroutine push_down(this, ind)
    type(pqueue_t), intent(inout) :: this
    integer, intent(in) :: ind
!
! Helper procedure to move items with lower priority down the heap.
!
    integer :: curi, chi, rchi

    curi = ind
    do while (curi < this%n/2 + 1)
      chi = 2*curi    ! select left child ...
      rchi = 2*curi+1
      if (rchi <= this%n) then
        ! ... or select the right child if it exists and has higher
        ! priority than the left child
        if (is_higher_priority(this%priorities(rchi),this%priorities(chi),this)) chi = rchi
      end if

      if (is_higher_priority(this%priorities(chi),this%priorities(curi),this)) then
        ! swap nodes if the selected child node has a higher priority than
        ! the current node
        call swap(this, chi, curi)
        curi = chi
      else
        exit
      end if
    end do
  end subroutine push_down


  function pqueue_insert(this, values, priority, ierr) result(handle)
    class(pqueue_t), intent(inout) :: this
    integer, intent(in) :: values(:), priority
    integer, intent(out), optional :: ierr
    type(handle_t) :: handle
!
! Insert (item, priority) pair to the queue. The "handle" is used to reference
! this item while it is in the queue. After item is removed from queue,
! the handle is no longer valid.
!
    if (size(values) /= this%chunksize) then
      call handle_error(ERR_INVALID_ARG_SIZE, 'pqueue_insert - unexpected size of newitem', ierr)
      return
    else if (.not. allocated(this%values)) then
      error stop 'pqueue_insert - pqueue is not itialized'
    else
      if (present(ierr)) ierr=ERR_OK
    end if

    call borrow_handle(this, handle)
    this%n = this%n + 1
    this%values(:,this%n) = values
    this%priorities(this%n) = priority
    this%handles(this%n) = handle
    this%hmap(handle%index_to_hmap) = this%n
    call bubble_up(this, this%n)
  end function pqueue_insert


  function pqueue_pop(this, top_priority, top_handle, ierr) result(top_value)
    class(pqueue_t), intent(inout) :: this
    integer, intent(out), optional :: top_priority
    type(handle_t), intent(out), optional :: top_handle
    integer, intent(out), optional :: ierr
    integer :: top_value(this%chunksize)
!
! Pop item with the highest priority (a lowest value of "p") from the queue.
! Throws an error if queue is empty.
! The item from the heap bottom is moved to the postion of removed item and
! then moved down the queue according its priority.
!
    if (.not. allocated(this%values)) then
      error stop 'pqueue_pop - pqueue is not itialized'
    else if (this%n < 1) then
      call handle_error(ERR_EMPTY, 'pop: empty queue', ierr)
      return
    else
      if (present(ierr)) ierr = ERR_OK
    end if

    ! swap top with the last item, then copy "top" into output variables
    if (this%n /= 1) call swap(this, 1, this%n)
    top_value = this%values(:, this%n)
    if (present(top_priority)) top_priority = this%priorities(this%n)
    if (present(top_handle)) top_handle = this%handles(this%n)

    ! unmark removed item from the hash-map
    call return_handle(this, this%handles(this%n))
    this%hmap(this%handles(this%n)%index_to_hmap) = HMAP_NULL

    ! push-down the last item (that is now the first element in the heap)
    this%n = this%n - 1
    if (this%n > 1) call push_down(this, 1)
  end function pqueue_pop


  subroutine pqueue_update_priority(this, handle, new_priority, ierr)
    class(pqueue_t), intent(inout) :: this
    type(handle_t), intent(in) :: handle
    integer, intent(in) :: new_priority
    integer, intent(out), optional :: ierr
!
! Update the priority of an item using handle.
!
    integer :: id, old_priority

    if (.not. allocated(this%values)) &
      & error stop 'pqueue_update - pqueue is not itialized'
    id = get_idheap(this, handle)
    if (id == HMAP_NULL) then
      call handle_error(ERR_INVALID_HANDLE, 'pqueue_update - invalid handle', ierr)
      return
    else
      if (present(ierr)) ierr = ERR_OK
    end if

    old_priority = this%priorities(id)
    this%priorities(id) = new_priority

    if (is_higher_priority(new_priority, old_priority, this)) then
      ! priority increased
      call bubble_up(this, id)
    else if (is_lower_priority(new_priority, old_priority, this)) then
      ! priority decreased
      call push_down(this, id)
    end if
  end subroutine pqueue_update_priority


  pure function pqueue_size(this) result(n)
    class(pqueue_t), intent(in) :: this
    integer :: n
    n = this%n
  end function pqueue_size


  pure function pqueue_empty(this) result(is_empty)
    class(pqueue_t), intent(in) :: this
    logical :: is_empty
    is_empty = this%n==0
  end function pqueue_empty



  ! -----
  ! OTHER
  ! -----

  subroutine handle_error(error_code, msg, ierr)
    integer, intent(in) :: error_code
    character(len=*), intent(in) :: msg
    integer, intent(inout), optional :: ierr

    write(efid,'("ERROR ",a)') msg
    if (present(ierr)) then
      ierr = error_code
    else
      error stop msg
    end if
  end subroutine handle_error

end module conts_mod
