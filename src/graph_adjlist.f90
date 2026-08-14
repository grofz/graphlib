! Copyright (C) 2026 Zdenek Grof
!
! This file is part of Graph library.
!
! Graph library is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published
! by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Graph library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with Graph library. If not, see <https://www.gnu.org/licenses/>.


! -----------------------------------------------------------------------------
! CONTAINS (only public objects listed)
!   type, public :: iterator_t
!   type, public :: adjlist_t
!   subroutine adjlist_initialize(this, capacity)
!   function adjlist_find(this, item) result(iterator)
!   function adjlist_contains(this, item)
!   subroutine adjlist_next(this, iterator, item)
!   function adjlist_has_next(this, iterator)
!   subroutine adjlist_add(this, item, skip_duplicity_check)
!   subroutine adjlist_remove(this, item, found_position)
! -----------------------------------------------------------------------------

  module graph_adjlist_mod
    implicit none (type, external)
    private

    type, public :: iterator_t
!     private
      integer :: i = 1
    end type
    interface iterator_t
      module procedure iterator_new
    end interface


    type, public :: adjlist_t
      private
      integer, allocatable :: arr(:)
      integer :: n = 0
    contains
      procedure :: initialize=>adjlist_initialize
      procedure :: find=>adjlist_find
      procedure :: contains=>adjlist_contains
      procedure :: next=>adjlist_next
      procedure :: next_noadvance=>adjlist_next_noadvance
      procedure :: advance=>adjlist_advance
      procedure :: has_next=>adjlist_has_next
      procedure :: back=>adjlist_back
      procedure :: add=>adjlist_add
      procedure :: remove=>adjlist_remove
      procedure :: size=>adjlist_size
    end type

    integer, parameter :: DEFAULT_CAPACITY = 5
    integer, parameter :: INCREASE_CAPACITY_MULTIPLIER = 2
    integer, parameter :: INCREASE_CAPACITY_ADDER = 0

    integer, parameter :: NULL_INDEX = huge(NULL_INDEX)

  contains

    pure function iterator_new() result(new)
      type(iterator_t) :: new
      new%i = 1
    end function iterator_new


    pure subroutine adjlist_initialize(this, capacity)
      class(adjlist_t), intent(inout) :: this
      integer, intent(in), optional :: capacity

      integer :: capacity0

      capacity0 = DEFAULT_CAPACITY
      if (present(capacity)) capacity0 = capacity

      ! allocate or reallocate "arr" to the given capacity
      if (allocated(this%arr)) then
        if (size(this%arr) /= capacity0) deallocate(this%arr)
      end if
      if (.not. allocated(this%arr)) allocate(this%arr(capacity0))

      this%n = 0
    end subroutine adjlist_initialize


    pure function adjlist_find(this, item) result(iterator)
      class(adjlist_t), intent(in) :: this
      integer, intent(in) :: item
      type(iterator_t) iterator

      integer :: i

      iterator%i = NULL_INDEX
      do i = 1, this%n
        if (this%arr(i)/=item) cycle
        iterator%i = i
        exit
      end do
    end function adjlist_find


    pure function adjlist_contains(this, item)
      class(adjlist_t), intent(in) :: this
      integer, intent(in) :: item
      logical adjlist_contains

      type(iterator_t) :: found

      found = this%find(item)
      adjlist_contains = found%i <= this%n
    end function adjlist_contains


    pure subroutine adjlist_next(this, iterator, item)
      class(adjlist_t), intent(in) :: this
      type(iterator_t), intent(inout) :: iterator
      integer, intent(out) :: item
!
! Return current item and advance iterator.
!
      if (allocated(this%arr)) then
        if (iterator%i <= this%n) then
          item = this%arr(iterator%i)
          iterator%i = iterator%i + 1
          return
        end if
      end if
      error stop 'adjlist_next - could not obtain next item'
    end subroutine adjlist_next


    pure subroutine adjlist_next_noadvance(this, iterator, item)
      class(adjlist_t), intent(in) :: this
      type(iterator_t), intent(in) :: iterator
      integer, intent(out) :: item
!
! Return current item but do not advance iterator.
!
      if (allocated(this%arr)) then
        if (iterator%i <= this%n) then
          item = this%arr(iterator%i)
          return
        end if
      end if
      error stop 'adjlist_next_noadvance - could not obtain next item'
    end subroutine adjlist_next_noadvance


    pure subroutine adjlist_advance(this, iterator)
      class(adjlist_t), intent(in) :: this
      type(iterator_t), intent(inout) :: iterator
!
! Advance iterator
!
      if (allocated(this%arr)) then
        if (iterator%i <= this%n) then
          iterator%i = iterator%i + 1
          return
        end if
      end if
      error stop 'adjlist_advance - error'
    end subroutine adjlist_advance


    pure function adjlist_has_next(this, iterator)
      class(adjlist_t), intent(in) :: this
      type(iterator_t), intent(in) :: iterator
      logical adjlist_has_next
!
! Return .true. if iterator points at the item and subroutine
! adjlist_next can be called.
!
      adjlist_has_next = iterator%i <= this%n
    end function adjlist_has_next


    pure subroutine adjlist_back(this, iterator, item)
      class(adjlist_t), intent(in) :: this
      type(iterator_t), intent(inout) :: iterator
      integer, intent(out), optional :: item
!
! Move iterator one item back and then return the item (if item provided)
! Should be called if the current item was deleted
!
      if (iterator%i > 1 .and. iterator%i <= this%n+1) then
        iterator%i = iterator%i - 1
        if (present(item)) item = this%arr(iterator%i)
        return
      else
        error stop 'adjlist_back - could not move one item back'
      end if
    end subroutine adjlist_back


    pure subroutine adjlist_add(this, item, skip_duplicity_check)
      class(adjlist_t), intent(inout) :: this
      integer, intent(in) :: item
      logical, intent(in), optional :: skip_duplicity_check

      logical :: skip_duplicity_check0

      skip_duplicity_check0 = .false.
      if (present(skip_duplicity_check)) skip_duplicity_check0 = skip_duplicity_check

      if (.not. allocated(this%arr)) &
          error stop 'adjlist_add - array not initialized'

      ! increase capacity if needed
      block
        integer, allocatable :: tmp(:)
        integer :: new_capacity
        if (this%n < size(this%arr)) then
          continue
        else if (this%n == size(this%arr)) then
          new_capacity = INCREASE_CAPACITY_MULTIPLIER*size(this%arr)+INCREASE_CAPACITY_ADDER
          allocate(tmp(new_capacity))
          tmp(1:this%n) = this%arr
          call move_alloc(tmp, this%arr)
        else
          error stop 'adjlist_add - something wrong with this%n and size of this%arr'
        end if
      end block

      ! avoid duplicit entries
      if (.not. skip_duplicity_check0) then
        if (this%contains(item)) &
            error stop 'adjlist_add - item is already present'
      end if

      this%n = this%n + 1
      this%arr(this%n) = item
    end subroutine adjlist_add


    pure subroutine adjlist_remove(this, item, found_position)
      class(adjlist_t), intent(inout) :: this
      integer, intent(in) :: item
      type(iterator_t), intent(in), optional :: found_position

      integer :: id
      type(iterator_t) :: found_position_here

      if (present(found_position)) then
        id = found_position%i
      else
        found_position_here = this%find(item)
        id = found_position_here%i
      end if

      ! verify "id" is sane
      if (id < 1 .or. id > this%n) &
          error stop 'adjlist_remove - invalid position given, or item not exists'
      if (this%arr(id) /= item) &
          error stop 'adjlist_remove - index points to an unexpected item'

      ! swap removed item with the last item in the array
      if (this%n > 1) this%arr(id) = this%arr(this%n)
      this%n = this%n - 1
    end subroutine adjlist_remove


    pure function adjlist_size(this) result(n)
      class(adjlist_t), intent(in) :: this
      integer n

      n = this%n
    end function adjlist_size

  end module graph_adjlist_mod
