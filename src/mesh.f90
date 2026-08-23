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
    use graph_adjlist_mod, only : adjlist_t, iterator_t
    use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR
    use conts_mod, only : queue_t, stack_t, INTEGER_MOLD
    implicit none (type, external)
    private

    ! Named local constants
    integer, parameter :: DEFAULT_CCAPACITY = 10, DEFAULT_PCAPACITY = 5

    ! Reference point used for 2D-grid only to define the positive orientation
    ! based on ordering the points.
    ! All point positions must be bellow the reference point.
    real(dp), parameter :: &
        ORIENTATION_2D_REFPOINT(3) = real([0.0, 0.0, 10.0], dp)

    ! Extending enumerators for "graph_handle_t"
    ! - type=1 for vertices and type=2 for edges in "src/graph.f90"
    integer(I1B), parameter :: POINT_HANDLE_TYPE = 3_I1B, &
        CELL_HANDLE_TYPE = 4_I1B, INVALID_HANDLE_TYPE= 0_I1B

    type, public :: point_t
      type(adjlist_t) :: depending_cells
      type(graph_handle_t) :: handle
      real(dp) :: position(3) = 0.0_dp
    end type point_t

    type, public :: cell_t
      type(graph_handle_t) :: points(4)
        ! ordered array of point handles
      type(graph_handle_t) :: ngbcells(4)
        ! array of cell handles to neighbouring cells
      type(graph_handle_t) :: dual_vertex
        ! handle to vertex in graph_t parent object
      type(graph_handle_t) :: handle
    contains
      procedure :: point_indices => cell_point_indices
    end type cell_t

    type, extends(graph_t), public :: mesh_t
      type(point_t), allocatable :: points(:)
      type(cell_t), allocatable :: cells(:)
      integer, allocatable, private :: pmap(:), cmap(:)
          ! storing position of points/cells in "points"/"cells" arrays
      integer :: npoints=NOT_INITIALIZED, ncells=NOT_INITIALIZED
      logical, private :: is_3d_mesh=.false.
          ! .false. = triangular mesh (cell defined from 3 points)
          ! .true. = tetrahedral mesh (cell defined from 4 points)
      type(queue_t), private :: free_phandles, free_chandles
    contains
      ! these procedures override procedures from graph_t class (note)
      procedure :: initialize => mesh_initialize
      procedure :: index_from_handle => mesh_index_from_handle
      procedure :: print => mesh_print
! TODO - override these or make them non-overridable in graph_t
     !procedure :: copy
     !procedure :: build_selection_masks
      procedure, non_overridable :: find_cell_id => mesh_find_cell_id
      procedure, non_overridable :: add_point => mesh_add_point
      procedure, non_overridable :: add_cell => mesh_add_cell
      procedure, non_overridable :: remove_cell => mesh_remove_cell
      procedure, non_overridable :: npoints_per_cell => mesh_npoints_per_cell
      procedure, non_overridable :: is_3d => mesh_is_3d
    end type mesh_t

  contains

    elemental function mesh_index_from_handle(this, handle) result(id)
      class(mesh_t), intent(in) :: this
      type(graph_handle_t), intent(in) :: handle
      integer id, index_to_map
!
! Return position of a point/cell in array using handle. If handle refers to
! the point/cell that is no longer in array, MAP_NULL is returned.
! If handle refers to vertex or edge, method of a parent class is used.
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
#ifdef DEBUG
        if (id /= MAP_NULL) then
          if (id<1 .or. id>this%npoints) error stop &
              'mesh_index_from_handle - point index out of bounds (internal error)'
        end if
#endif
      case(CELL_HANDLE_TYPE)
        if (index_to_map > 0 .and. index_to_map <= size(this%cmap)) then
          id = this%cmap(index_to_map)
          if (id/=MAP_NULL) then
            ! verify version matches the stored one
            if (.not. (this%cells(id)%handle==handle)) id = MAP_NULL
          end if
        end if
#ifdef DEBUG
        if (id /= MAP_NULL) then
          if (id<1 .or. id>this%ncells) error stop &
              'mesh_index_from_handle - cell index out of bounds (internal error)'
        end if
#endif
      case default
        ! delegate to parent class
        id = this%graph_t%index_from_handle(handle)
      end select
    end function mesh_index_from_handle


    subroutine borrow_mesh_handle(this, handle_type, handle)
      class(mesh_t), intent(inout) :: this
      integer(i1b), intent(in) :: handle_type
      type(graph_handle_t), intent(out) :: handle

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
        error stop &
            'borrow_mesh_handle: invalid handle type, expecting point/cell type only'
      end select
    end subroutine borrow_mesh_handle


    subroutine return_mesh_handle(this, handle)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: handle

      type(graph_handle_t) :: reused_handle

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

    pure integer function mesh_npoints_per_cell(this) result(n)
      class(mesh_t), intent(in) :: this
      n = 3
      if (this%is_3d_mesh) n = 4
    end function mesh_npoints_per_cell


    pure logical function mesh_is_3d(this) result(is)
      class(mesh_t), intent(in) :: this
      is = this%is_3d_mesh
    end function mesh_is_3D


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

      ! initialize graph_t components
      call this%graph_t%initialize(vcapacity=vcapacity, ecapacity=ecapacity, &
          is_directed_graph=.false.)

      ! two-dimensional or three dimensional?
      this%is_3d_mesh = .false. ! 2d is a default at the moment
      if (present(is_3d)) this%is_3d_mesh = is_3d

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
        type(graph_handle_t) :: handle
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
          'increase_points_capacity - new_capacity <= old_capacity '

      ! reallocate "points" and "pmap"
      allocate(tmp_points(new_capacity0))
      allocate(tmp_map(new_capacity0), source=MAP_NULL)
      tmp_points(1:old_capacity) = this%points
      tmp_map(1:old_capacity) = this%pmap
      call move_alloc(tmp_points, this%points)
      call move_alloc(tmp_map, this%pmap)

      block ! create fresh handles
        integer :: i
        type(graph_handle_t) :: new_handle
        do i=old_capacity+1, new_capacity0
          new_handle = graph_handle_t(id=i, version=1, type=POINT_HANDLE_TYPE)
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
          'increase_cells_capacity - new_capacity <= old_capacity '

      ! reallocate "cells" and "cmap"
      allocate(tmp_cells(new_capacity0))
      allocate(tmp_map(new_capacity0), source=MAP_NULL)
      tmp_cells(1:old_capacity) = this%cells
      tmp_map(1:old_capacity) = this%cmap
      call move_alloc(tmp_cells, this%cells)
      call move_alloc(tmp_map, this%cmap)

      block ! create fresh handles
        integer :: i
        type(graph_handle_t) :: new_handle
        do i=old_capacity+1, new_capacity0
          new_handle = graph_handle_t(id=i, version=1, type=CELL_HANDLE_TYPE)
          call this%free_chandles%enqueue(transfer(new_handle, INTEGER_MOLD))
        end do
      end block
    end subroutine increase_cells_capacity


    function mesh_add_point(this, position) result(handle)
      class(mesh_t), intent(inout) :: this
      real(dp), intent(in) :: position(3)
      type(graph_handle_t) :: handle
!
! Add point to the mesh.
!
      if (.not. this%is_initialized()) then
        error stop 'mesh_add_point - not initialized'
      end if

      ! For a 2-D mesh, reject positions above the reference point
      if (.not. this%is_3d_mesh) then
        if (position(3) >= ORIENTATION_2D_REFPOINT(3)) error stop &
            'mesh_add_point - expecting z-coordinate bellow reference point'
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


    subroutine relocate_point(this, handle, newid)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: handle
      integer, intent(in) :: newid

      integer :: oldid

      if (handle%get_handle_type() /= POINT_HANDLE_TYPE) &
          error stop 'relocate_point - wrong handle type'
      oldid = this%index_from_handle(handle)
      if (oldid == MAP_NULL) &
          error stop 'relocate_point - point no more exists'
      if (newid < 1 .or. newid > this%npoints) &
          error stop 'relocate_point - newid out of bounds'

      ! copy point and update record in "pmap"
      this%points(newid) = this%points(oldid)
      this%pmap(handle%get_index_to_map()) = newid
    end subroutine relocate_point


    function mesh_add_cell(this, point_handles) result(handle)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: point_handles(4)
      type(graph_handle_t) :: handle
!
! Add mesh cell. Dual vertex and edges connecting neighbouring vertices are
! also added.
!
      integer :: n, i, j, pids(4), ngb_pids(4)
      integer, allocatable :: found_cids(:)

      n = this%npoints_per_cell()

      if (.not. this%is_initialized()) then
        error stop 'mesh_add_cell - not initialized'
      else if (any(point_handles(1:n)%get_handle_type() /= POINT_HANDLE_TYPE)) then
        error stop 'mesh_add_cell - point type handles expected'
      end if

      ! All points must exist and be unique
      block
        integer :: pids0(4)
        logical :: unique

        pids0(1:n) = this%index_from_handle(point_handles(1:n))
        if (any(pids0(1:n)==MAP_NULL)) &
            error stop 'mesh_add_cell - a point not present (invalid handle)'
        unique = .true. ! innocent until found guilty
        do i = 1, n-1
          do j = i+1, n
            if (pids0(i)==pids0(j)) unique = .false.
          end do
        end do
        if (.not. unique) error stop 'mesh_add_cell - all points must be unique'
      end block

      ! Order points
      pids = order_point_indices(this, point_handles)

      ! Check if cell already exists
      found_cids = this%find_cell_id(pids(1:n))
      select case(size(found_cids))
      case(0)
        continue ! all ok, can proceed
      case(1)
        error stop 'mesh_add_cell - cell already exists'
      case default
        error stop 'mesh_add_cell - too many cells (internal error)'
      end select
      deallocate(found_cids)

      ! We verified there is no cell among points
      ! Set components of added cell
      call borrow_mesh_handle(this, CELL_HANDLE_TYPE, handle)
      this%ncells = this%ncells + 1
      associate (new_cell => this%cells(this%ncells))
        new_cell%handle = handle

        do i = 1, n
          new_cell%points(i) = this%points(pids(i))%handle

          ! Set connection with neighbouring cells
          ! A neighbouring cell across point i must contain all points except
          ! point i
          found_cids = this%find_cell_id([pids(1:i-1), pids(i+1:n)])
          select case (size(found_cids))
          case(0)
            ! no oposite cell accros point "i" / explicitly set to null
            new_cell%ngbcells(i) = &
              graph_handle_t(id=MAP_NULL, type=CELL_HANDLE_TYPE, version=1)
print '("Point ",i0," - no ngb across")', this%index_from_handle(new_cell%points(i))
          case(1)
            ! set connection
            associate(ngb_cell => this%cells(found_cids(1)))
              new_cell%ngbcells(i) = ngb_cell%handle

              ! which node in neigbouring cell is accross the added cell?
              ngb_pids = ngb_cell%point_indices(this)
              do j=1, n
                if (all(pids(1:n)/=ngb_pids(j))) exit
              end do
              if (j==n+1) error stop &
                  'mesh_add_cell - opposite point in ngb cell not found (internal error)'
              ! verify that handle, we are about to set, points to null
              if (ngb_cell%ngbcells(j)%get_index_to_map()/=MAP_NULL) &
                  error stop 'mesh_add_cell - a connection exists (internal error)'
              ! set the connection back
              ngb_cell%ngbcells(j) = handle
print '("Point ",i0," - across is cell ",i0)', this%index_from_handle(new_cell%points(i)), this%index_from_handle(ngb_cell%handle)
            end associate
          case default
            error stop 'mesh_add_cell - more than one neighbouring cell found (internal error)'
          end select
          deallocate(found_cids)
        end do

        ! Make sure handles at unused position are set to null in 2D-meshes
        if (.not. this%is_3d_mesh) then
          new_cell%points(4) = &
              graph_handle_t(id=MAP_NULL, type=POINT_HANDLE_TYPE, version=1)
          new_cell%ngbcells(4) = &
              graph_handle_t(id=MAP_NULL, type=CELL_HANDLE_TYPE, version=1)
        end if

        ! Add a dual vertex and, if appropriate, also edges representing
        ! connection between neighboring cells
        block
          type(graph_handle_t) :: edge
          integer :: ngb_id, v_ipar(VSIZE_IPAR), e_ipar(ESIZE_IPAR)
          real(dp) :: v_rpar(VSIZE_RPAR), e_rpar(ESIZE_RPAR)

          new_cell%dual_vertex = this%add_vertex(v_ipar, v_rpar)
          do i=1, n
            ngb_id = this%index_from_handle(new_cell%ngbcells(i))
            if (ngb_id==MAP_NULL) cycle
            edge = this%add_edge(new_cell%dual_vertex, &
                this%cells(ngb_id)%dual_vertex, e_ipar, e_rpar)
          end do
        end block

      end associate
      this%cmap(handle%get_index_to_map()) = this%ncells

      ! Add the new cell to the points's list...
      do i = 1, n
        call this%points(pids(i))%depending_cells%add(this%ncells)
      end do
print '("Cell added. Handle is ",i0)', this%index_from_handle(handle)

    end function mesh_add_cell


    subroutine mesh_remove_cell(this, handle)
      class(mesh_t), intent(inout) :: this
      class(graph_handle_t), intent(in) :: handle
!
! TODO Documentation
!
      integer :: icell, pids(4), i, n

      if (handle%get_handle_type()/=CELL_HANDLE_TYPE) error stop &
          'mesh_remove_cell - invalid handle type, cell type expected'

      icell = this%index_from_handle(handle)
      if (icell == MAP_NULL) error stop &
          'mesh_remove_cell - cell no longer exists'

      n = this%npoints_per_cell()

      ! Remove cell refrence from list(s) of its points.
      ! It is ok if a point no longer exist, but if it exists, the
      ! reference to the cell must be present in its list.
      pids = this%cells(icell)%point_indices(this)
      do i=1, n
        if (pids(i) /= MAP_NULL) &
            call this%points(pids(i))%depending_cells%remove(icell)
      end do

      ! Remove dual-vertex (and associated edges)
      call this%remove_vertex(this%cells(icell)%dual_vertex)

      ! Unlink removed cell from neighbouring cells
      block
        integer :: ngbid, ngb_pids(4), j

        do i=1, this%npoints_per_cell()
          ngbid = this%index_from_handle(this%cells(icell)%ngbcells(i))
          if (ngbid==MAP_NULL) cycle
          associate(ngb_cell => this%cells(ngbid))
            ! which node in neigbouring cell is accross the removed cell?
            ngb_pids = ngb_cell%point_indices(this)
            do j=1, n
              if (all(pids(1:n)/=ngb_pids(j))) exit
            end do
            if (j==n+1) error stop &
                'mesh_remove_cell - opposite point in ngb cell not found (internal error)'
            ! verify that handle, we are about to clear, points to icell
            if (this%index_from_handle(ngb_cell%ngbcells(j))/=icell) &
                error stop 'mesh_remove_cell - unexpected connection link (internal error)'
            ! unlink the connection back
            ngb_cell%ngbcells(j) = &
                graph_handle_t(id=MAP_NULL, type=CELL_HANDLE_TYPE, version=1)
          end associate
        end do
      end block

      ! Nullify cmap entry and reuse the handle
      this%cmap(handle%get_index_to_map()) = MAP_NULL
      call return_mesh_handle(this, handle)

      ! Relocate the last cell to the "hole" after removed cell 
      if (icell /= this%ncells) then
        call relocate_cell(this, this%cells(this%ncells)%handle, icell)
      end if
      this%ncells = this%ncells - 1
      
    end subroutine mesh_remove_cell


    subroutine relocate_cell(this, handle, newid)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: handle
      integer, intent(in) :: newid

      integer :: oldid, i, pids(4)

      if (handle%get_handle_type() /= CELL_HANDLE_TYPE) &
          error stop 'relocate_cell - wrong handle type'
      oldid = this%index_from_handle(handle)
      if (oldid == MAP_NULL) &
          error stop 'relocate_cell - cell no more exists'
      if (newid < 1 .or. newid > this%ncells) &
          error stop 'relocate_cell - newid out of bounds'

      ! copy cell and update record in "cmap"
      this%cells(newid) = this%cells(oldid)
      this%cmap(handle%get_index_to_map()) = newid

      ! update depending cells lists of respective points
      pids = this%cells(newid)%point_indices(this)
      do i=1, this%npoints_per_cell()
        if (pids(i)/=MAP_NULL) call update_list(this%points(pids(i))%depending_cells)
      end do

    contains
      subroutine update_list(list) ! internal procedure
        type(adjlist_t), intent(inout) :: list

        type(iterator_t) :: found_oldid
        if (list%contains(newid)) &
            error stop 'relocate cell - newid present in list would lead to duplicity'
        found_oldid = list%find(oldid)
        if (.not. list%has_next(found_oldid)) &
            error stop 'relocate cell - oldid not found in list'
        call list%remove(oldid, found_oldid)
        call list%add(newid, skip_duplicity_check=.true.)
      end subroutine
    end subroutine relocate_cell


    subroutine mesh_print(this, fid)
      class(mesh_t), intent(in) :: this
      integer, intent(in) :: fid
!
! Print all mesh data for debugging
!
      character(len=:), allocatable :: text_meshdim
      integer :: i, j, cid
      type(iterator_t) :: iterator

      if (.not. this%is_initialized()) then
        write(fid,'("Mesh not initialized (WARNING)")')
        return
      end if

      if (this%is_3d_mesh) then
        text_meshdim = '3D mesh'
      else
        text_meshdim = '2D mesh'
      end if
      write(fid,'("--- mesh dump ---")')
      write(fid, '(a," with ",i0," points and ",i0," cells")') &
          text_meshdim, this%npoints, this%ncells

      ! information about points
      do i=1, this%npoints
        write(fid, '("P-",i0," position",3(1x,g0))') &
          this%points(i)%handle%get_index_to_map(), this%points(i)%position
      end do
      write(fid, *)
      do i=1, this%npoints
        write(fid, '("P-",i0," incident cells")',advance='no') &
            this%points(i)%handle%get_index_to_map()
        iterator = iterator_t()
        do while (this%points(i)%depending_cells%has_next(iterator))
          call this%points(i)%depending_cells%next(iterator, cid)
          write(fid, '(1x,i0)', advance='no') &
            this%cells(cid)%handle%get_index_to_map()
        end do
        write(fid,*)
      end do

      ! information about cells
      do i=1, this%ncells
        write(fid,'(/,"C-",i0,/,"--from")',advance='no') &
          this%cells(i)%handle%get_index_to_map()
        do j=1, this%npoints_per_cell()
          write(fid,'(1x,"P-",i0)',advance='no') &
            this%cells(i)%points(j)%get_index_to_map()
        end do
        write(fid,'(/,"--connected with cells")',advance='no')
        do j=1, this%npoints_per_cell()
          write(fid,'(1x,i0)',advance='no') &
            this%cells(i)%ngbcells(j)%get_index_to_map()
        end do
        write(fid,'(/,"--dual with V-",i0)') &
          this%cells(i)%dual_vertex%get_index_to_map()
      end do

      call this%graph_t%print(fid)
      write(fid,'("--- end of mesh dump ---",/)')

    end subroutine mesh_print


    ! ---------------------
    ! Mesh helper functions
    ! ---------------------


!TODO when pure stack%pop is ready, make this function pure
    function mesh_find_cell_id(this, pids) result(cids)
      class(mesh_t), intent(in) :: this
      integer, intent(in) :: pids(:)
      integer, allocatable :: cids(:)
!
! Given a group of 1 to 3 (or 4) points, find all cells associated with these
! points and return an array of cell indices. Return empty sized array if no
! such cell exists.
!
      type(iterator_t) :: iterator
      type(stack_t) :: found_cells
      integer :: cid_current, pids_current(4), i, j
      logical :: has_point(4)

      ! Validate input
      if (size(pids) > this%npoints_per_cell()) error stop &
          'mesh_find_cell_id - too many points given'
      if (any(pids<1 .or. pids>this%npoints)) error stop &
          'mesh_find_cell_id - point indices out of bounds'
      call found_cells%initialize(chunksize=size(transfer(1,INTEGER_MOLD)),capacity=5)

      ! it should be enough to look in depcells list of a single point only
      if (size(pids)>0) then
        iterator = iterator_t()
        do while (this%points(pids(1))%depending_cells%has_next(iterator))
          call this%points(pids(1))%depending_cells%next(iterator, cid_current)
          pids_current = this%cells(cid_current)%point_indices(this)

          ! does current cell consists of all given points?
          has_point = .false.
          has_point(1) = .true.
          do i = 2, size(pids)
            do j=1, this%npoints_per_cell()
              if (pids(i)==pids_current(j)) then
                has_point(i) = .true.
                exit
              end if
            end do
            if (.not. has_point(i)) exit
          end do
          if (all(has_point(1:size(pids)))) &
              call found_cells%push(transfer(cid_current,INTEGER_MOLD))
        end do ! next cell from the list
      end if

      ! prepare output array
      allocate(cids(found_cells%size()))
      i = 0
      do while(.not. found_cells%empty())
        i = i + 1
        cids(i) = transfer(found_cells%pop(), 1)
      end do
      if (i /= size(cids)) error stop &
          'mesh_find_cell_id - stack consumation irregularity (internal error)'
    end function mesh_find_cell_id


   !pure function mesh_find_mirror_cell_id(this, cell, ploc) result(cid)
   !  class(mesh_t), intent(in) :: this
   !  type(cell_t), intent(in) :: cell
   !end function mesh_find_mirror_cell_id


    ! ------------
    ! Cell methods
    ! ------------
    pure function cell_point_indices(this, mesh) result(ids)
      class(cell_t), intent(in) :: this
      type(mesh_t), intent(in) :: mesh
      integer :: ids(4)

      integer :: i

      if (.not. mesh%is_3d_mesh) ids(4) = MAP_NULL
      do i=1, mesh%npoints_per_cell()
        ids(i) = mesh%index_from_handle(this%points(i))
      end do
    end function cell_point_indices


    pure function order_point_indices(this, points) result(pids)
      class(mesh_t), intent(in) :: this
      type(graph_handle_t), intent(in) :: points(4)
      integer :: pids(4)
!
! Order points for positive orientation.
! 
      real(dp), parameter :: eps = 10 * epsilon(1.0_dp)
      real(dp), parameter :: p_ref(3) = ORIENTATION_2D_REFPOINT
      real(dp) :: d, tol
      integer :: n, itmp
      type(graph_handle_t) :: points0(4), ptmp

      n = this%npoints_per_cell()

      ! The first two positions will be points with the lowest index_to_map.
      points0 = points
      associate (loc=>minloc(points0(1:n)%get_index_to_map(), dim=1))
        if (loc /= 1) then
          ptmp = points0(1)
          points0(1) = points0(loc)
          points0(loc) = ptmp
        end if
      end associate
      associate (loc=>minloc(points0(2:n)%get_index_to_map(), dim=1))
        if (loc /= 1) then
          ptmp = points0(2)
          points0(2) = points0(loc+1)
          points0(loc+1) = ptmp
        end if
      end associate

      ! Point indices in the actual mesh
      pids(1:n) = this%index_from_handle(points0(1:n))
      if (any(pids(1:n)<1 .or. any(pids(1:n)>this%npoints))) error stop &
          'order_point_indices - point indices out of bounds (internal error)'

      ! Orientation: det [(p2-p1); (p3-p1); (p4-p1)] > 0
      ! For 2d mesh, an arbitrary reference point is used instead of p4
      block
        real(dp) :: a(3), b(3), axb(3), c(3)

        ! 2D-mesh: although not used, p4 must be associated with valid item
        if (.not. this%is_3d_mesh) pids(4) = pids(1)
        associate( &
          p1=>this%points(pids(1))%position, &
          p2=>this%points(pids(2))%position, &
          p3=>this%points(pids(3))%position, &
          p4=>this%points(pids(4))%position)

          a = p2 - p1
          b = p3 - p1
          if (this%is_3d_mesh) then
            c = p4 - p1
          else
            c = p_ref - p1
          end if
        end associate
        ! axb is a normal vector to the p1-p2-p3 plane
        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)
        ! d determines the side of the plane on which p4/p_ref lies
        d = dot_product(axb, c)
        ! tolerance for the signed volume calculation
        tol = eps*max(1.0_dp, maxval(abs(a))*maxval(abs(b))*maxval(abs(c)))
      end block

      if (abs(d)<tol) then
        error stop 'order_point_indices - degenerate positions'
      elseif (d < -tol) then
        itmp = pids(2)
        pids(2) = pids(3)
        pids(3) = itmp
      else
        continue
      end if
    end function order_point_indices

  end module mesh_mod
