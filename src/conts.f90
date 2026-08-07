! =============================================================================
! Container data structures
!
! General-purpose container data structures used by the library.
! The containers store integer arrays internally, so arbitrary data can
! be stored by using transfer() as a wrapper. The array size must be set
! during initialization (see "chunksize" argument) and stays until container
! is reinitialized.
!
! These containers are provided
!
! - stack_t
! LIFO stack with dynamically growing storage.
!
! - queue_t
! FIFO queue with dynamically growing storage.
!
! - pqueue_t
! Priority queue implemented as an indexed binary heap. Items are
! identified by handles, allowing their priority to be queried, updated,
! or removed in O(1) / O(log N) time as appropriate.
!
! The priority queue maintains a handle for each stored item. Handles remain
! valid while the corresponding item is in the queue and become invalid after
! the item is removed. This allows algorithms such as Dijkstra's shortest-path
! search to update priorities without searching the queue for the item.
! =============================================================================

module conts_mod
  use, intrinsic :: iso_fortran_env, only : dp=>real64
  implicit none (type, external)
  private

  integer, parameter :: NOT_INITIALIZED = -1
  integer, parameter :: DEFAULT_CAPACITY = 10
  integer, parameter :: INTEGER_MOLD(0) = [integer ::]

  type, abstract :: container_t
    integer :: n = NOT_INITIALIZED
    integer, allocatable :: values(:,:)
  contains
    procedure, non_overridable :: size => container_size
    procedure, non_overridable :: empty => container_empty
    procedure, non_overridable :: initialized => container_initialized
    procedure(container_clear_ai), deferred :: clear
  end type

  abstract interface
    pure subroutine container_clear_ai(this)
      import container_t
      implicit none
      class(container_t), intent(inout) :: this
    end subroutine
  end interface


  ! STACK
  type, extends(container_t), public :: stack_t
    private
  contains
    procedure :: initialize => stack_initialize
    procedure :: push => stack_push
    procedure :: pop => stack_pop
    procedure :: export => stack_export
    procedure :: peek => stack_peek
    procedure :: clear => stack_clear
  end type


  ! QUEUE
  type, extends(container_t), public :: queue_t
    private
    integer :: rear
  contains
    procedure :: initialize => queue_initialize
    procedure :: enqueue => queue_enqueue
    procedure :: dequeue => queue_dequeue
    procedure :: export => queue_export
    procedure :: peek => queue_peek
    procedure :: clear => queue_clear
  end type


  ! PRIORITY QUEUE
  integer, parameter :: HMAP_NULL = -1
  integer, parameter, public :: &
    PQUEUE_MIN = 1, & ! lower P is higher priority
    PQUEUE_MAX = 2    ! higher P is higher priority

  type, public :: handle_t
    private
    integer :: index_to_hmap = HMAP_NULL
    integer :: version = 1
  contains
    procedure, private :: handle_eq, handle_write_formatted
    generic :: operator(==) => handle_eq
    ! There may be a compiler error: with this on, the allocatable array
    ! of handle_t does not deallocate automatically and must be deallocated
    ! manually. We keep it commented out for now as this is not needed
    !generic :: write(formatted) => handle_write_formatted
  end type
  interface handle_t
    module procedure handle_new
  end interface

  type, extends(container_t), public :: pqueue_t
    private
    ! "values" ... heap of values inherited from base class
    real(DP), allocatable :: priorities(:)
      ! heap of priorities
    type(handle_t), allocatable :: handles(:)
      ! heap of handles
    integer, allocatable :: hmap(:)
      ! handle to map to item's position in heaps
    type(queue_t) :: free_handles
    integer :: ordering = PQUEUE_MIN
  contains
    procedure :: initialize => pqueue_initialize
    procedure :: insert => pqueue_insert
    procedure :: pop => pqueue_pop
    procedure :: remove => pqueue_remove
    procedure :: priority => pqueue_priority
    procedure :: update_priority => pqueue_update_priority
    procedure :: contains => pqueue_contains
    procedure :: peek => pqueue_peek
    procedure :: export => pqueue_export
    procedure :: export_priorities => pqueue_export_priorities
    procedure :: export_handles => pqueue_export_handles
    procedure :: valid => pqueue_valid
    procedure :: clear => pqueue_clear
  end type
  ! TODO: heapify, update (insert/update_priority)

contains

! --------------
! COMMON METHODS
! --------------

  pure function container_size(this) result(n)
    class(container_t), intent(in) :: this
    integer :: n
    n = this%n
  end function container_size


  pure function container_empty(this) result(is_empty)
    class(container_t), intent(in) :: this
    logical :: is_empty
    is_empty = this%n < 1
  end function container_empty


  pure function container_initialized(this) result(is)
    class(container_t), intent(in) :: this
    logical :: is
    is = this%n /= NOT_INITIALIZED
  end function container_initialized


  pure subroutine stack_clear(this)
    class(stack_t), intent(inout) :: this
    if (.not. this%initialized()) error stop 'clear - container not initialized'
    this%n = 0
  end subroutine stack_clear


  pure subroutine queue_clear(this)
    class(queue_t), intent(inout) :: this
    if (.not. this%initialized()) error stop 'clear - container not initialized'
    this%n = 0
    this%rear = 1
  end subroutine queue_clear


  pure subroutine pqueue_clear(this)
    class(pqueue_t), intent(inout) :: this
    integer :: i
    if (.not. this%initialized()) error stop 'clear - container not initialized'
    do i=1, this%n
      this%hmap(this%handles(i)%index_to_hmap) = HMAP_NULL
      call return_handle(this, this%handles(i))
    end do
    this%n = 0
  end subroutine pqueue_clear


! -------------
! STACK / QUEUE
! -------------

  pure subroutine stack_initialize(this, chunksize, capacity)
    class(stack_t), intent(inout) :: this
    integer, intent(in), optional :: chunksize, capacity

    integer :: new_capacity, new_chunksize

    if (allocated(this%values)) deallocate(this%values)
    if (present(chunksize)) then
      new_chunksize = chunksize
    else
      new_chunksize = 1
    end if
    if (present(capacity)) then
      new_capacity = max(1, capacity)
    else
      new_capacity = DEFAULT_CAPACITY
    end if
    this%n = 0
    allocate(this%values(new_chunksize,0))
    call stack_increase_capacity(this, new_capacity)
  end subroutine


  pure subroutine queue_initialize(this, chunksize, capacity)
    class(queue_t), intent(inout) :: this
    integer, intent(in), optional :: chunksize, capacity

    integer :: new_capacity, new_chunksize

    if (allocated(this%values)) deallocate(this%values)
    if (present(chunksize)) then
      new_chunksize = chunksize
    else
      new_chunksize = 1
    end if
    if (present(capacity)) then
      new_capacity = max(1, capacity)
    else
      new_capacity = DEFAULT_CAPACITY
    end if
    this%n = 0
    this%rear = 1
    allocate(this%values(new_chunksize,0))
    call queue_increase_capacity(this, new_capacity)
  end subroutine


  pure subroutine stack_increase_capacity(this, new_capacity)
    class(stack_t), intent(inout) :: this
    integer, intent(in) :: new_capacity

    integer, allocatable :: tmp(:,:)

    if (new_capacity <= this%n) &
      & error stop 'stack_increase_capacity - internal error'
    allocate(tmp(size(this%values,dim=1), new_capacity))
    tmp(:,1:this%n) = this%values
    call move_alloc(tmp, this%values)
  end subroutine stack_increase_capacity


  pure subroutine queue_increase_capacity(this, new_capacity)
    class(queue_t), intent(inout) :: this
    integer, intent(in) :: new_capacity

    integer, allocatable :: tmp(:,:)
    integer :: front, c1, c2, old_capacity

    if (new_capacity <= this%n) &
      & error stop 'queue_increase_capacity - internal error'
    old_capacity = size(this%values,dim=2)
    allocate(tmp(size(this%values,dim=1), new_capacity))
    if (old_capacity /= 0) then
      front = modulo(this%rear-this%n-1, old_capacity) + 1
    else
      front = 1
    end if
    if (front >= this%rear .and. this%n/=0) then
      c1 = old_capacity-front+1
      c2 = this%n-c1
      tmp(:,1:c1) = this%values(:,front:old_capacity)
      tmp(:,c1+1:c1+c2) = this%values(:,1:c2)
    else
      tmp(:,1:this%n) = this%values(:,front:this%rear-1)
    end if
    call move_alloc(tmp, this%values)
    this%rear = this%n+1
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


  pure subroutine stack_push(this, newitem)
    class(stack_t), intent(inout) :: this
    integer, intent(in) :: newitem(:)

    if (.not. this%initialized()) then
      error stop 'stack_push - stack is not initialized'
    else if (size(newitem) /= size(this%values,dim=1)) then
      error stop 'stack_push - newitem size invalid'
    end if
    if (size(this%values,dim=2) == this%n) &
      & call stack_increase_capacity(this, 2*size(this%values,dim=2))

    ! push new item
    this%n = this%n + 1
    this%values(:,this%n) = newitem
  end subroutine stack_push


  pure subroutine queue_enqueue(this, newitem)
    class(queue_t), intent(inout) :: this
    integer, intent(in) :: newitem(:)

    if (.not. this%initialized()) then
      error stop 'queue_enqueue - queue is not itialized'
    else if (size(newitem) /= size(this%values,dim=1)) then
      error stop 'queue_enqueue - newitem size invalid'
    end if
    if (size(this%values,dim=2) == this%n) &
      & call queue_increase_capacity(this, 2*size(this%values,dim=2))

    ! enqueue new item
    this%values(:,this%rear) = newitem
    this%n = this%n + 1
    this%rear = modulo(this%rear, size(this%values,dim=2))+1
  end subroutine queue_enqueue


  function stack_pop(this) result(pop_item)
    class(stack_t), intent(inout) :: this
    integer :: pop_item(size(this%values,dim=1))

    if (this%n < 1) error stop 'stack_pop - empty stack'
    pop_item = this%values(:,this%n)
    this%n = this%n - 1
  end function stack_pop


  function queue_dequeue(this) result(pop_item)
    class(queue_t), intent(inout) :: this
    integer :: pop_item(size(this%values,dim=1))

    if (this%n < 1) error stop 'queue_dequeue - empty queue'
    pop_item = this%values(:, modulo(this%rear-this%n-1, size(this%values,dim=2)) + 1)
    this%n = this%n - 1
  end function queue_dequeue


  function stack_peek(this) result(item)
    class(stack_t), intent(in) :: this
    integer :: item(size(this%values,dim=1))

    if (this%n < 1) error stop 'stack_peek - empty stack'
    item = this%values(:,this%n)
  end function stack_peek


  function queue_peek(this) result(item)
    class(queue_t), intent(in) :: this
    integer :: item(size(this%values,dim=1))

    if (this%n < 1) error stop 'queue_peek - empty queue'
    item = this%values(:, modulo(this%rear-this%n-1, size(this%values,dim=2)) + 1)
  end function queue_peek


  pure function stack_export(this) result(values)
    class(stack_t), intent(in) :: this
    integer :: values(size(this%values,dim=1), this%n)
    values = this%values(:,1:this%n)
  end function


  pure function queue_export(this) result(values)
    class(queue_t), intent(in) :: this
    integer :: values(size(this%values,dim=1), this%n)

    integer :: front, c1, c2

    if (size(this%values,dim=2) /= 0) then
      front = modulo(this%rear-this%n-1, size(this%values,dim=2)) + 1
    else
      front = 1
    end if
    if (front >= this%rear .and. this%n/=0) then
      c1 = size(this%values,dim=2)-front+1
      c2 = this%n-c1
      values(:,1:c1) = this%values(:,front:size(this%values,dim=2))
      values(:,c1+1:c1+c2) = this%values(:,1:c2)
    else
      values(:,1:this%n) = this%values(:,front:this%rear-1)
    end if
  end function


  ! ------------------------------
  ! PRIORITY QUEUE implementation
  ! ------------------------------

  pure function handle_eq(a, b) result(eq)
    class(handle_t), intent(in) :: a, b
    logical :: eq
    eq = a%version==b%version .and. a%index_to_hmap==b%index_to_hmap
  end function handle_eq


  pure function handle_new() result(new)
    type(handle_t) :: new
    new%version = 1
    new%index_to_hmap = HMAP_NULL
  end function handle_new


  pure function is_higher_priority(pa, pb, this) result(is)
    real(dp), intent(in) :: pa, pb
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
  end function is_higher_priority


  pure function is_lower_priority(pa, pb, this) result(is)
    real(dp), intent(in) :: pa, pb
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
  end function is_lower_priority


  pure subroutine pqueue_initialize(this, chunksize, capacity, ordering)
    class(pqueue_t), intent(inout) :: this
    integer, intent(in), optional :: chunksize, capacity, ordering

    integer :: new_capacity, new_chunksize
    type(handle_t) :: handle

    if (present(chunksize)) then
      new_chunksize = chunksize
    else
      new_chunksize = 1
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
    allocate(this%values(new_chunksize,0))
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


  pure subroutine pqueue_increase_capacity(this, new_capacity)
    class(pqueue_t), intent(inout) :: this
    integer, intent(in), optional :: new_capacity

    integer :: old_capacity, new_capacity0
    integer, allocatable :: tmp2(:,:), tmp1(:)
    real(dp), allocatable :: tmp3(:)
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

    allocate(tmp3(new_capacity0))
    tmp3(1:old_capacity) = this%priorities
    call move_alloc(tmp3, this%priorities)

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


  pure subroutine return_handle(this, handle)
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
    real(dp) :: tmp_real

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
    tmp_real = this%priorities(i)
    this%priorities(i) = this%priorities(j)
    this%priorities(j) = tmp_real

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


  function pqueue_insert(this, values, priority) result(handle)
    class(pqueue_t), intent(inout) :: this
    integer, intent(in) :: values(:)
    real(dp), intent(in) :: priority
    type(handle_t) :: handle
!
! Insert (item, priority) pair to the queue. The "handle" is used to reference
! this item while it is in the queue. After item is removed from queue,
! the handle is no longer valid.
!
    if (.not. this%initialized()) then
      error stop 'pqueue_insert - pqueue is not itialized'
    else if (size(values) /= size(this%values,dim=1)) then
      error stop 'pqueue_insert - invalid size of newitem'
    end if

    call borrow_handle(this, handle)
    this%n = this%n + 1
    this%values(:,this%n) = values
    this%priorities(this%n) = priority
    this%handles(this%n) = handle
    this%hmap(handle%index_to_hmap) = this%n
    call bubble_up(this, this%n)
  end function pqueue_insert


  function pqueue_pop(this, top_priority, top_handle) result(top_value)
    class(pqueue_t), intent(inout) :: this
    real(dp), intent(out), optional :: top_priority
    type(handle_t), intent(out), optional :: top_handle
    integer :: top_value(size(this%values,dim=1))
!
! Pop item with the highest priority (a lowest value of "p") from the queue.
! Throws an error if queue is empty.
! The item from the heap bottom is moved to the postion of removed item and
! then moved down the queue according its priority.
!
    if (this%n < 1) error stop 'pop: empty pqueue'

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


  pure subroutine pqueue_remove(this, handle)
    class(pqueue_t), intent(inout) :: this
    type(handle_t), intent(in) :: handle
!
! Remove an item from the queue
!
    integer :: id

    if (.not. this%initialized()) error stop &
        'pqueue_remove - pqueue is not itialized'
    id = get_idheap(this, handle)
    if (id == HMAP_NULL) error stop 'pqueue_remove - invalid handle'

    ! swap removed item with the last item
    if (this%n /= id) call swap(this, id, this%n)
    ! unmark removed item from the hash-map
    call return_handle(this, this%handles(this%n))
    this%hmap(this%handles(this%n)%index_to_hmap) = HMAP_NULL
    ! push-down the last item (that took place of removed element in the heap)
    this%n = this%n - 1
    if (this%n > 1) call push_down(this, id)
  end subroutine pqueue_remove


  function pqueue_peek(this, top_priority, top_handle) result(top_value)
    class(pqueue_t), intent(in) :: this
    real(dp), intent(out), optional :: top_priority
    type(handle_t), intent(out), optional :: top_handle
    integer :: top_value(size(this%values,dim=1))
!
! Look at item with highest priority in the queue without removing it
!
    if (this%n < 1) error stop 'peek: empty pqueue'

    ! swap top with the last item, then copy "top" into output variables
    top_value = this%values(:, 1)
    if (present(top_priority)) top_priority = this%priorities(1)
    if (present(top_handle)) top_handle = this%handles(1)
  end function pqueue_peek


  pure function pqueue_priority(this, handle) result(priority)
    class(pqueue_t), intent(in) :: this
    type(handle_t), intent(in) :: handle
    real(dp) :: priority
!
! Return priority of the item
!
    integer :: id

    if (.not. this%initialized()) error stop 'pqueue_priority - not initialized'
    id = get_idheap(this, handle)
    if (id == HMAP_NULL) error stop  'pqueue_priority - invalid handle'
    priority = this%priorities(id)

  end function pqueue_priority


  subroutine pqueue_update_priority(this, handle, new_priority)
    class(pqueue_t), intent(inout) :: this
    type(handle_t), intent(in) :: handle
    real(dp), intent(in) :: new_priority
!
! Update the priority of an item using handle.
!
    integer :: id
    real(dp) :: old_priority

    if (.not. this%initialized()) error stop 'pqueue_update - not itialized'
    id = get_idheap(this, handle)
    if (id == HMAP_NULL) error stop 'pqueue_update - invalid handle'

    old_priority = this%priorities(id)
    this%priorities(id) = new_priority

    if (is_higher_priority(new_priority, old_priority, this)) then
      ! priority increased
      call bubble_up(this, id)
    else if (is_lower_priority(new_priority, old_priority, this)) then
      ! priority decreased
      call push_down(this, id)
    else
      error stop 'pqueue_update_priority - priority seems unchanged'
    end if
  end subroutine pqueue_update_priority


  pure function pqueue_export(this) result(values)
    class(pqueue_t), intent(in) :: this
    integer :: values(size(this%values,dim=1), this%n)
    values = this%values(:,1:this%n)
  end function pqueue_export


  pure function pqueue_export_priorities(this) result(priorities)
    class(pqueue_t), intent(in) :: this
    real(dp) :: priorities(this%n)
    priorities = this%priorities(1:this%n)
  end function pqueue_export_priorities


  pure function pqueue_export_handles(this) result(handles)
    class(pqueue_t), intent(in) :: this
    type(handle_t) :: handles(this%n)
    handles = this%handles(1:this%n)
  end function pqueue_export_handles


  pure recursive function validate_heap(this, ind) result(valid)
    class(pqueue_t), intent(in) :: this
    integer, intent(in) :: ind
    logical :: valid

    logical :: left_valid, right_valid

    left_valid = .true.
    right_valid = .true.
    if (2*ind <= this%n) left_valid = validate_heap(this, 2*ind) .and. .not. &
      & is_higher_priority(this%priorities(2*ind), this%priorities(ind), this)
    if (2*ind+1 <= this%n) right_valid = validate_heap(this, 2*ind+1) .and. .not. &
      & is_higher_priority(this%priorities(2*ind+1), this%priorities(ind), this)
    valid = left_valid .and. right_valid
  end function validate_heap


  pure function pqueue_valid(this) result(valid)
    class(pqueue_t), intent(in) :: this
    logical :: valid

    if (.not. this%initialized()) then
      ! uninitialized quueue is assumed valid
      valid = .true.
      return
    end if

    ! capacity = size + free_handles
    valid = size(this%values, dim=2) == this%n + this%free_handles%size()

    ! active handles
    block
      integer :: i
      do i=1, this%n
        if (this%hmap(this%handles(i)%index_to_hmap)/=i) then
          valid = .false.
          exit
        end if
      enddo
    end block

    ! null items in hmap
    valid = valid .and. count(this%hmap==HMAP_NULL) == this%free_handles%size()

    ! priority heap is correct
    valid = valid .and. validate_heap(this, 1)
  end function pqueue_valid


  ! -----
  ! OTHER
  ! -----

  subroutine handle_write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)
    class(handle_t), intent(in) :: dtv
    integer, intent(in) :: unit
    character(len=*), intent(in) :: iotype
    integer, intent(in) :: v_list(:)
    integer, intent(out) :: iostat
    character(len=*), intent(inout) :: iomsg

    ! just to avoig "unused dummy variables warning'
    block
      logical :: dummy
      dummy = iotype==''
      if (size(v_list)>0) dummy=v_list(1)==1
    end block

    iostat = 0
    iomsg = ''
    write(unit, '(I0," (v",I0,")")') dtv%index_to_hmap, dtv%version
  end subroutine handle_write_formatted

end module conts_mod
