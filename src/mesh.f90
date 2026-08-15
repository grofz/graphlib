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


  module mesh_mod
    use iso_fortran_env, only : dp=>real64, I1B=>int8
    use graph_mod, only : graph_t, graph_handle_t=>handle_t, MAP_NULL, &
        NOT_INITIALIZED
    use graph_adjlist_mod, only : adjlist_t
    use conts_mod, only : queue_t
    implicit none (type, external)
    private

    ! Named local constants
    integer, parameter :: DEFAULT_CCAPACITY = 10, DEFAULT_PCAPACITY = 5
    integer, parameter :: INTEGER_MOLD(0) = [integer ::]

    ! type=1 for vertices and type=2 for edges in "src/graph.f90"
    integer(I1B), parameter :: POINT_HANDLE_TYPE = 3_I1B, &
        CELL_HANDLE_TYPE = 4_I1B, INVALID_HANDLE_TYPE= 0_I1B


    type, extends(graph_handle_t), public :: mesh_handle_t
    contains
      procedure :: get_index_to_map => mhandle_get_index_to_map
    end type mesh_handle_t

    type, public :: point_t
      type(adjlist_t) :: depending_cells
      type(mesh_handle_t) :: handle
      real(dp) :: position(3) = 0.0_dp
    end type point_t

    type :: cell_t
      type(mesh_handle_t) :: points(4)
      type(mesh_handle_t) :: ngbcells(4)
      type(mesh_handle_t) :: handle
    contains
      procedure :: point_indices => cell_point_indices
    end type cell_t

    type, extends(graph_t), public :: mesh_t
      type(point_t), allocatable :: points(:)
      type(cell_t), allocatable :: cells(:)
      integer, allocatable, private :: pmap(:), cmap(:)
          ! storing position of points/cells in "points"/"cells" arrays
      integer :: npoints=NOT_INITIALIZED, ncells=NOT_INITIALIZED
      logical, private :: is_3d=.false.
          ! .false. = triangular mesh (cell defined from 3 points)
          ! .true. = tetrahedral mesh (cell defined from 4 points)
      type(queue_t), private :: free_phandles, free_chandles
    contains
      ! these procedures override procedures from graph_t class (note)
      procedure :: initialize => mesh_initialize
      procedure :: get_index_from_handle => mesh_get_index_from_handle
! TODO - override these or make them non-overridable in graph_t
     !procedure :: copy
     !procedure :: build_selection_masks
     !procedure :: print
      procedure, non_overridable :: add_point => mesh_add_point
    end type mesh_t

  contains

    pure function mesh_get_index_from_handle(this, handle) result(id)
      class(mesh_t), intent(in) :: this
      type(graph_handle_t), intent(in) :: handle
      integer id, index_to_map
!
! Return position of a point/cell in array using handle. If handle refers to
! the point/cell that is no longer in array, MAP_NULL is returned.
! If handle refers to vertex or edge, method of a parent class is used
!
      id = MAP_NULL
      index_to_map = handle%get_index_to_map()
      select case(handle%get_handle_type())
      case(POINT_HANDLE_TYPE)
        if (index_to_map > 0 .and. index_to_map <= size(this%pmap)) then
          id = this%pmap(index_to_map)
          if (id/=MAP_NULL) then
            ! verify version matches the stored one
            if (.not. (this%points(id)%handle==handle)) id = MAP_NULL
          end if
        end if
      case(CELL_HANDLE_TYPE)
        if (index_to_map > 0 .and. index_to_map <= size(this%cmap)) then
          id = this%cmap(index_to_map)
          if (id/=MAP_NULL) then
            ! verify version matches the stored one
            if (.not. (this%cells(id)%handle==handle)) id = MAP_NULL
          end if
        end if
      case default
        ! delegate to parent class
        id = this%graph_t%get_index_from_handle(handle)
      end select
    end function mesh_get_index_from_handle


    elemental integer function mhandle_get_index_to_map(this, graph) result(id)
      class(mesh_handle_t), intent(in) :: this
      class(graph_t), intent(in), optional :: graph
      if (present(graph)) then
        select type(graph)
        class is (mesh_t)
          id = mesh_get_index_from_handle(graph, this%graph_handle_t)
        class default
          error stop 'mhandle_get_index_to_map - mesh_handle_t object must be used with mesh_t object only'
        end select
      else
        ! delegate to the getter from the parent class
        id = this%graph_handle_t%get_index_to_map()
      end if
    end function mhandle_get_index_to_map


    subroutine borrow_mesh_handle(this, handle_type, handle)
      class(mesh_t), intent(inout) :: this
      integer(i1b), intent(in) :: handle_type
      type(mesh_handle_t), intent(out) :: handle

      select case(handle_type)
      case(POINT_HANDLE_TYPE)
        if (this%free_phandles%size()==0) call increase_points_capacity(this)
        if (this%free_phandles%size()==0) error stop 'borrow_mesh_handle - no more P-handles available'
        handle = transfer(this%free_phandles%dequeue(), handle)
      case(CELL_HANDLE_TYPE)
        if (this%free_chandles%size()==0) call increase_cells_capacity(this)
        if (this%free_chandles%size()==0) error stop 'borrow_mesh_handle - no more C-handles available'
        handle = transfer(this%free_chandles%dequeue(), handle)
      case default
        error stop 'borrow_mesh_handle: invalid handle_type'
      end select
    end subroutine borrow_mesh_handle


    subroutine return_mesh_handle(this, handle)
      class(mesh_t), intent(inout) :: this
      type(mesh_handle_t), intent(in) :: handle

      type(mesh_handle_t) :: reused_handle

      reused_handle = handle
      call reused_handle%advance_version()
      select case(reused_handle%get_handle_type())
      case(POINT_HANDLE_TYPE)
        call this%free_phandles%enqueue(transfer(reused_handle,INTEGER_MOLD))
      case(CELL_HANDLE_TYPE)
        call this%free_chandles%enqueue(transfer(reused_handle,INTEGER_MOLD))
      case default
        error stop 'return_mesh_handle: handle provided has invalid handle_type'
      end select
    end subroutine return_mesh_handle


    ! ----------------------
    ! Mesh basic operations
    ! ----------------------
    subroutine mesh_initialize(this, vcapacity, ecapacity, is_directed_graph, &
        pcapacity, ccapacity, is_3d)
      class(mesh_t), intent(inout) :: this
      integer, intent(in), optional :: pcapacity, ccapacity, vcapacity, ecapacity
      logical, intent(in), optional :: is_3d, is_directed_graph

      ! graph is always indirected for a mesh_t object
      if (present(is_directed_graph)) then
        if (is_directed_graph) error stop &
            'mesh_initialize - mesh_t object requires undirected graph'
      end if

      ! two-dimensional or three dimensional?
      this%is_3d = .false. ! 2d is a default at the moment
      if (present(is_3d)) this%is_3d = is_3d

      ! reallocate all arrays to zero size
      if (allocated(this%points)) deallocate(this%points)
      allocate(this%points(0))
      if (allocated(this%cells)) deallocate(this%cells)
      allocate(this%cells(0))
      if (allocated(this%pmap)) deallocate(this%pmap)
      allocate(this%pmap(0))
      if (allocated(this%cmap)) deallocate(this%cmap)
      allocate(this%cmap(0))

      this%npoints = 0
      this%ncells = 0

      block ! initialize queues of free handles
        type(mesh_handle_t) :: handle
        call this%free_phandles%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
        call this%free_chandles%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
      end block

      block ! allocate initial capacity
        integer :: new_capacity
        new_capacity = DEFAULT_PCAPACITY
        if (present(pcapacity)) new_capacity = pcapacity
        call increase_points_capacity(this, new_capacity)
        new_capacity = DEFAULT_CCAPACITY
        if (present(ccapacity)) new_capacity = ccapacity
        call increase_cells_capacity(this, new_capacity)
      end block

      ! initialize graph_t components
      call this%graph_t%initialize(vcapacity=vcapacity, ecapacity=ecapacity, &
          is_directed_graph=.false.)
    end subroutine mesh_initialize


    subroutine increase_points_capacity(this, new_capacity)
      class(mesh_t), intent(inout) :: this
      integer, intent(in), optional :: new_capacity

      integer :: old_capacity, new_capacity0
      type(point_t), allocatable :: tmp_points(:)
      integer, allocatable :: tmp_map(:)

      old_capacity = size(this%points)
      if (present(new_capacity)) then
        new_capacity0 = new_capacity
      else
        new_capacity0 = 2*old_capacity
      end if

      if (new_capacity0 <= old_capacity) error stop &
          'increase_points_capacity - new_capaciry <= old_capacity '

      ! reallocate "points" and "pmap"
      allocate(tmp_points(new_capacity0))
      allocate(tmp_map(new_capacity0), source=MAP_NULL)
      tmp_points(1:old_capacity) = this%points
      tmp_map(1:old_capacity) = this%pmap
      call move_alloc(tmp_points, this%points)
      call move_alloc(tmp_map, this%pmap)

      block ! create fresh handles
        integer :: i
        type(mesh_handle_t) :: new_handle
        do i=old_capacity+1, new_capacity0
          new_handle%graph_handle_t = &
              graph_handle_t(id=i, version=1, type=POINT_HANDLE_TYPE)
          call this%free_phandles%enqueue(transfer(new_handle, INTEGER_MOLD))
        end do
      end block
    end subroutine increase_points_capacity


    subroutine increase_cells_capacity(this, new_capacity)
      class(mesh_t), intent(inout) :: this
      integer, intent(in), optional :: new_capacity

      integer :: old_capacity, new_capacity0
      type(cell_t), allocatable :: tmp_cells(:)
      integer, allocatable :: tmp_map(:)

      old_capacity = size(this%cells)
      if (present(new_capacity)) then
        new_capacity0 = new_capacity
      else
        new_capacity0 = 2*old_capacity
      end if

      if (new_capacity0 <= old_capacity) error stop &
          'increase_cells_capacity - new_capaciry <= old_capacity '

      ! reallocate "edges" and "vmap"
      allocate(tmp_cells(new_capacity0))
      allocate(tmp_map(new_capacity0), source=MAP_NULL)
      tmp_cells(1:old_capacity) = this%cells
      tmp_map(1:old_capacity) = this%cmap
      call move_alloc(tmp_cells, this%cells)
      call move_alloc(tmp_map, this%cmap)

      block ! create fresh handles
        integer :: i
        type(mesh_handle_t) :: new_handle
        do i=old_capacity+1, new_capacity0
          new_handle%graph_handle_t = &
              graph_handle_t(id=i, version=1, type=CELL_HANDLE_TYPE)
          call this%free_chandles%enqueue(transfer(new_handle, INTEGER_MOLD))
        end do
      end block
    end subroutine increase_cells_capacity


    function mesh_add_point(this, position) result(handle)
      class(mesh_t), intent(inout) :: this
      real(dp), intent(in) :: position(3)
      type(mesh_handle_t) :: handle

      if (.not. this%is_initialized()) then
        error stop 'mesh_add_point - not initialized'
      end if

      call borrow_mesh_handle(this, POINT_HANDLE_TYPE, handle)
      this%npoints = this%npoints + 1
      associate (new_point => this%points(this%npoints))
        new_point%position = position
        call new_point%depending_cells%initialize()
        new_point%handle = handle
      end associate
      this%pmap(handle%get_index_to_map()) = this%npoints
    end function mesh_add_point


    ! -------------------
    ! Cell TPB procedures
    ! -------------------
    pure function cell_point_indices(this, mesh) result(ids)
      class(cell_t), intent(in) :: this
      type(mesh_t), intent(in) :: mesh
      integer :: ids(4)

      integer :: i
      do i=1,4
        ids(i) = mesh%get_index_from_handle(this%ngbcells(i)%graph_handle_t)
      end do
    end function cell_point_indices

  end module mesh_mod