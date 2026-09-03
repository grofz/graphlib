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
        NOT_INITIALIZED, is_vertex_selected, is_edge_selected, conjugate_gradient
    use graph_adjlist_mod, only : adjlist_t, iterator_t
    use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR
    use conts_mod, only : queue_t, stack_t, INTEGER_MOLD
    implicit none (type, external)
    private
public integrate_pde, solve_3x3 ! TODO for testing temporarily

    ! Named local constants
    integer, parameter :: DEFAULT_CCAPACITY = 10, DEFAULT_PCAPACITY = 5
    real(dp), parameter :: GEOMETRY_TOL_FACTOR = 1000.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: NUMERICAL_TOL = 10.0_dp*epsilon(1.0_dp)

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
      type(graph_handle_t) :: ngb_cells(4)
        ! array of cell handles to neighbouring cells
      type(graph_handle_t) :: dual_vertex
        ! handle to vertex in graph_t parent object
      type(graph_handle_t) :: handle
    contains
      procedure :: point_indices => cell_point_indices
      procedure :: geometry => cell_geometry
    end type cell_t

    type, public :: cell_geometry_t
      real(dp) :: centre(3) ! circumcentre
      real(dp) :: volume
      real(dp) :: area_vector(3,4) ! outward face/edge area vectors
      real(dp) :: face_distance(4) ! centre-to-face distance
    end type cell_geometry_t

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
      procedure :: npoints_per_cell => mesh_npoints_per_cell
! TODO - override these or make them non-overridable in graph_t
     !procedure :: copy
     !procedure :: build_selection_masks
      procedure, non_overridable :: find_cell_id => mesh_find_cell_id
      procedure, non_overridable :: add_point => mesh_add_point
      procedure, non_overridable :: add_cell => mesh_add_cell
      procedure, non_overridable :: remove_point => mesh_remove_point
      procedure, non_overridable :: remove_cell => mesh_remove_cell
      procedure, non_overridable :: is_3d => mesh_is_3d
      procedure, non_overridable :: append_rectilinear_mesh => mesh_append_rectilinear_mesh
    end type mesh_t

  contains

    elemental function mesh_index_from_handle(this, handle) result(id)
      class(mesh_t), intent(in) :: this
      type(graph_handle_t), intent(in) :: handle
      integer id, index_to_map
!
! Given an object handle, return a valid[*] array position of the object,
! or MAP_NULL if handle refers to an object no longer present in the array.
! Points and cells are processed here; handles refering to vertices or
! edges are processed by the parent class method graph_index_from_handle().
!
! [*]  Post-call validation that the returned index is within the active
!      object array is unnecessary.
!
      id = MAP_NULL
      index_to_map = handle%get_index_to_map()

      select case(handle%get_handle_type())
      case(POINT_HANDLE_TYPE)
        if (index_to_map > 0 .and. index_to_map <= size(this%pmap)) then
          id = this%pmap(index_to_map)
          if (id/=MAP_NULL) then
            ! Verify that the given handle matches the stored one.
            if (id < 1 .or. id > this%npoints) error stop &
                'mesh_index_from_handle - point index out of bounds (internal error)'
            if (.not. (this%points(id)%handle==handle)) id = MAP_NULL
          end if
        end if

      case(CELL_HANDLE_TYPE)
        if (index_to_map > 0 .and. index_to_map <= size(this%cmap)) then
          id = this%cmap(index_to_map)
          if (id/=MAP_NULL) then
            ! Verify that the given handle matches the stored one.
            if (id < 1 .or. id > this%ncells) error stop &
                'mesh_index_from_handle - cell index out of bounds (internal error)'
            if (.not. (this%cells(id)%handle==handle)) id = MAP_NULL
          end if
        end if

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
    end function mesh_is_3d


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
          new_handle = graph_handle_t( &
              id=i, version=this%get_graph_id(), type=POINT_HANDLE_TYPE)
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
          new_handle = graph_handle_t( &
              id=i, version=this%get_graph_id(), type=CELL_HANDLE_TYPE)
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


    subroutine mesh_remove_point(this, handle)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: handle
!
! Remove point. Cells associated with the point will be also removed.
!
      integer :: ipoint

      if (handle%get_handle_type() /= POINT_HANDLE_TYPE) error stop &
          'mesh_remove_point - invalid handle type, point type handle expected'
      ipoint = this%index_from_handle(handle)
      if (ipoint == MAP_NULL) &
          error stop 'mesh_remove_point - point no longer present in graph'

      ! Automatically remove all associated cells.
      block
        type(iterator_t) :: iterator
        integer :: icell
        do
          iterator = iterator_t()
          if (.not. this%points(ipoint)%depending_cells%has_next(iterator)) exit
          call this%points(ipoint)%depending_cells%next(iterator, icell)
          call mesh_remove_cell(this, this%cells(icell)%handle)
        end do
      end block

      ! Defensive
#ifdef DEBUG
      if (this%points(ipoint)%depending_cells%size()>0) &
          error stop 'mesh_remove_point - could not remove associated cells'
#endif

      ! Nullify pmap and return handle
      this%pmap(handle%get_index_to_map()) = MAP_NULL
      call return_mesh_handle(this, handle)

      ! Relocate the last point to fill "hole" after removed point
      if (ipoint /= this%npoints) then
        call relocate_point(this, this%points(this%npoints)%handle, ipoint)
      end if
      this%npoints = this%npoints - 1
    end subroutine mesh_remove_point


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


    function mesh_add_cell(this, point_handles, position_id) result(handle)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: point_handles(4)
      integer, intent(in), optional :: position_id
      type(graph_handle_t) :: handle
!
! Add mesh cell. Dual vertex and edges connecting neighbouring vertices are
! also added.
!
      integer :: n, i, j, pids(4), ngb_pids(4)
      integer, allocatable :: found_cids(:)

      if (present(position_id)) then
        if (position_id < 1 .or. position_id+2 > VSIZE_RPAR) error stop &
          'mesh_add_cell - position_id out of bounds'
      end if

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
            new_cell%ngb_cells(i) = &
              graph_handle_t(id=MAP_NULL, type=CELL_HANDLE_TYPE, version=1)
          case(1)
            ! set connection
            associate(ngb_cell => this%cells(found_cids(1)))
              new_cell%ngb_cells(i) = ngb_cell%handle

              ! which node in neigbouring cell is accross the added cell?
              ngb_pids = ngb_cell%point_indices(this)
              do j=1, n
                if (all(pids(1:n)/=ngb_pids(j))) exit
              end do
              if (j==n+1) error stop &
                  'mesh_add_cell - opposite point in ngb cell not found (internal error)'
              ! verify that handle, we are about to set, points to null
              if (ngb_cell%ngb_cells(j)%get_index_to_map()/=MAP_NULL) &
                  error stop 'mesh_add_cell - a connection exists (internal error)'
              ! set the connection back
              ngb_cell%ngb_cells(j) = handle
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
          new_cell%ngb_cells(4) = &
              graph_handle_t(id=MAP_NULL, type=CELL_HANDLE_TYPE, version=1)
        end if

        ! Add a dual vertex and, if appropriate, also edges representing
        ! connection between neighboring cells
        block
          type(graph_handle_t) :: edge
          integer :: ngb_id, v_ipar(VSIZE_IPAR), e_ipar(ESIZE_IPAR)
          real(dp) :: v_rpar(VSIZE_RPAR), e_rpar(ESIZE_RPAR)

          v_rpar = 0.0_dp
          e_rpar = 0.0_dp
          v_ipar = 0
          e_ipar = 0
          if (present(position_id)) then
            ! position of dual vertex (stored for visualization)
            associate (x=>v_rpar(position_id:position_id+2))
              if (this%is_3d()) then
                x = circumcentre_sphere( &
                    this%points(pids(1))%position, &
                    this%points(pids(2))%position, &
                    this%points(pids(3))%position, &
                    this%points(pids(4))%position )
              else
                x = circumcentre_circle( &
                    this%points(pids(1))%position, &
                    this%points(pids(2))%position, &
                    this%points(pids(3))%position )
              end if
            end associate
          end if
          new_cell%dual_vertex = this%add_vertex(v_ipar, v_rpar)
          do i=1, n
            ngb_id = this%index_from_handle(new_cell%ngb_cells(i))
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

    end function mesh_add_cell


    subroutine mesh_remove_cell(this, handle)
      class(mesh_t), intent(inout) :: this
      type(graph_handle_t), intent(in) :: handle
!
! Remove cell from the mesh. The underlying graph that represents cells
! connectivity is updated accordingly and the vertex dual to the cell is
! also removed.
!
      integer :: icell, pids(4), i, n

      if (handle%get_handle_type()/=CELL_HANDLE_TYPE) error stop &
          'mesh_remove_cell - invalid handle type, cell type expected'

      icell = this%index_from_handle(handle)
      if (icell == MAP_NULL) error stop &
          'mesh_remove_cell - cell no longer exists'

      n = this%npoints_per_cell()

      ! Remove cell reference from list(s) of its points.
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
          ngbid = this%index_from_handle(this%cells(icell)%ngb_cells(i))
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
            if (this%index_from_handle(ngb_cell%ngb_cells(j))/=icell) &
                error stop 'mesh_remove_cell - unexpected connection link (internal error)'
            ! unlink the connection back
            ngb_cell%ngb_cells(j) = &
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
            this%cells(i)%ngb_cells(j)%get_index_to_map()
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
    pure function cell_point_indices(this, mesh, null_allowed) result(ids)
      class(cell_t), intent(in) :: this
      type(mesh_t), intent(in) :: mesh
      logical, intent(in), optional :: null_allowed
      integer :: ids(4)
!
! Return (and optionally validate) indices of points forming the cell.
!
! For a 2D-mesh, ids(4) = ids(1).
!
! Unless "null_allowed = .true." is given, all returned indices are checked to
! be valid (and no post-function validation is required).
!
      integer :: i
      logical :: null_allowed0

      do i=1, mesh%npoints_per_cell()
        ids(i) = mesh%index_from_handle(this%points(i))
      end do

      ! By default, all points are expected to exist
      null_allowed0 = .false.
      if (present(null_allowed)) null_allowed0 = null_allowed
      if (.not. null_allowed0) then
        if (any(ids(1:mesh%npoints_per_cell())<1 ) .or. &
            any(ids(1:mesh%npoints_per_cell())>mesh%npoints)) error stop &
            'cell_point_indices - a point not in mesh'
      end if

      ! For a 2D mesh, ids(4) is usually not used
      ! Just to have a valid index for all array items
      if (.not. mesh%is_3d_mesh) ids(4) = ids(1)
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
        ! (based on mean term contributing to the determinant)
        tol = eps/6 * &
            (abs(a(1)*b(2)*c(3)) + abs(a(1)*b(3)*c(2)) + &
             abs(a(2)*b(1)*c(3)) + abs(a(2)*b(3)*c(1)) + &
             abs(a(3)*b(1)*c(2)) + abs(a(3)*b(2)*c(1)))
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


    function cell_geometry(this, mesh) result(geom)
!   pure function cell_geometry(this, mesh) result(geom)
      class(cell_t), intent(in) :: this
      type(mesh_t), intent(in) :: mesh
      type(cell_geometry_t) :: geom
!
! TODO Documentation
!
      logical :: is_3d
      integer :: pids(4), n
      real(dp) :: avec_123(3), tol
      real(dp), target :: positions(3,4)
      real(dp), pointer :: p1(:), p2(:), p3(:), p4(:)

      ! Get point positions
      is_3D = mesh%is_3d()
      n = mesh%npoints_per_cell()
      pids = this%point_indices(mesh)
      positions(:,1) = mesh%points(pids(1))%position
      positions(:,2) = mesh%points(pids(2))%position
      positions(:,3) = mesh%points(pids(3))%position
      positions(:,4) = mesh%points(pids(4))%position
      p1 => positions(:,1)
      p2 => positions(:,2)
      p3 => positions(:,3)
      p4 => positions(:,4) ! not used in a 2D mesh

      ! Assert non-degenerate positions
      tol = maxval(abs(p2-p1))
      tol = max(tol, maxval(abs(p3-p1)))
      tol = max(tol, maxval(abs(p4-p1)))
      tol = tol * GEOMETRY_TOL_FACTOR
      select case(n)
      case(3) ! a 2D mesh
        if (points_colinear(p1, p2, p3, tol)) error stop &
          'vell_geometry - 2D cell is co-linear'
      case(4) ! a 3D mesh
        if (points_coplanar(p1, p2, p3, p4, tol)) error stop &
          'cell_geometry - 3D cell is co-planar'
      case default
        error stop 'cell_geometry - invalid npoints_per_cell'
      end select

      ! Area vector of triangle base
      avec_123 = triangle_area_vector(p1, p2, p3)

      ! Cell volume
      if (is_3d) then
        geom%volume = dot_product(p4-p1, avec_123) / 3.0_dp
#ifdef DEBUG
        if (geom%volume <= 0.0_dp) error stop &
          'cell_geometry - negative volume'
#endif
      else
        geom%volume = norm2(avec_123)
      end if

      ! Cell centre
      if (is_3d) then
        geom%centre = circumcentre_sphere(p1, p2, p3, p4)
      else
        geom%centre = circumcentre_circle(p1, p2, p3)
      end if

      ! Face area vectors / Face distance
      block
        integer :: iface
        real(dp), target :: selected_pos(3,3)
        real(dp) :: fopp(3)
        real(dp), pointer :: f1(:), f2(:), f3(:)

        do iface = 1, n
          selected_pos = face_positions(positions, iface, n)
          f1 => selected_pos(:,1)
          f2 => selected_pos(:,2)
          f3 => selected_pos(:,3) ! unused for a 2D grid
          fopp = positions(:,iface)

          associate(avec => geom%area_vector(:,iface), &
              fdist => geom%face_distance(iface))
            if (is_3d) then
              avec = triangle_area_vector(f1, f2, f3)
              ! Area vector now points in the positive direction of the plane
              ! defined by f1-f2-f3 orientation. Area vector points outwards
              ! the cell, if the vertex opposite to the face is located on the
              ! negative side of the plane.
              ! Flip area vector direction if needed.
              if (point_plane_distance(f1,f2,f3,fopp) > 0.0_dp) avec = -avec

              fdist = abs(point_plane_distance(f1, f2, f3, geom%centre))
            else
              avec = vec_product(f2-f1, avec_123) / norm2(avec_123)

              fdist = point_line_distance(f1, f2, geom%centre)
            end if
          end associate
        end do
      end block

      ! Sum of area vectors should be zero
#ifdef DEBUG
      block
        real(dp), parameter :: SUM_TOL_FACTOR = 100.0_dp*epsilon(1.0_dp)
        real(dp) :: avec_sum(3), scale
        integer :: i
        avec_sum = sum(geom%area_vector(:,1:n), dim=2)
        scale = 0.0
        do i = 1, n
          scale = max(scale, norm2(geom%area_vector(:,i)))
        end do
        if (norm2(avec_sum) >= SUM_TOL_FACTOR * scale) error stop &
          'cell_geometry - area vectors do not sum up to zero'
      end block
#endif

    end function cell_geometry


    pure function face_positions(original_x, iface, npoints_per_cell) result(x)
      real(dp), intent(in) :: original_x(3,4)
      integer, intent(in) :: iface
      integer, intent(in) :: npoints_per_cell
      real(dp) :: x(3,3)
!
! Helper function - return "x" with "original_x(:,iface)" removed,
!
!    2D mesh:                  3D mesh:
!    iface  output             iface  output
!        1  x(:,2), x(:,3)         1  x(:,2), x(:,3), x(:,4)
!        2  x(:,3), x(:,1)         2  x(:,3), x(:,4), x(:,1)
!        3  x(:,1), x(:,2)         3  x(:,4), x(:,1), x(:,2)
!                                  4  x(:,1), x(:,2), x(:,3)
      integer :: i, orders(4)

      if (npoints_per_cell /= 3 .and. npoints_per_cell /=4) &
          error stop 'face_positions - npoints_per_cell 3 or 4 expected'
      if (iface < 1 .or. iface > npoints_per_cell) &
          error stop 'face_positions - iface is invalid'

      orders = [0, 1, 2, 3]
      orders = mod(orders+(iface-1), npoints_per_cell)

      do i = 2, npoints_per_cell
        x(:,i-1) = original_x(:,orders(i)+1)
      end do
    end function face_positions


    ! -------------
    ! Generate mesh  TODO temporarily here
    ! -------------

    subroutine mesh_append_rectilinear_mesh(this, p0, p1, p2, cell_size, &
        rel_shift, position_id, boundary, boundary_offset)
      class(mesh_t), intent(inout) :: this
      real(dp), intent(in) :: p0(3), p1(3), p2(3)
      real(dp), intent(in) :: cell_size, rel_shift
      integer, intent(in), optional :: position_id
      type(graph_handle_t), allocatable, intent(out) :: boundary(:)
      integer, intent(out) :: boundary_offset(0:4)
!
! TODO Documentation
!
! IN
!   this        - Mesh object
!   p0, p1, p2  - Corners of box defining the mesh border
!   cell_size   - Maximum length of the triangle side.
!   rel_shift   - Randomized pertrubation of grid points
!   position_id - (optional) To store position of dual and ghost vertices in the
!                 underlying vertex-edge graph. Can be used for visualisation
!                 and debugging.
!
! OUT
!   boundary - Handles to ghost-vertices connected to the boundary cells of
!              the mesh.
!   boundary_offset - identifies which handles belong to which mesh side:
!     side 01 - boundary( boundary_offset(0)+1 : boundary_offset(1) )
!     side 02 - boundary( boundary_offset(1)+1 : boundary_offset(2) )
!     side 23 - boundary( boundary_offset(2)+1 : boundary_offset(3) )
!     side 13 - boundary( boundary_offset(3)+1 : boundary_offset(4) )
!
!               (3)
!     P2 +---------------+ P3
!        |               |
!    (2) |               | (4)
!        |               |
!     P0 +---------------+ P1
!               (1)
!
      real(dp) :: base1(3), base2(3), tile_size(2)
      integer :: i, j
      type(graph_handle_t), allocatable :: p_corners(:,:), p_mids(:,:)
      type(queue_t) :: bitems_01, bitems_02, bitems_23, bitems_13
      type(graph_handle_t) :: p(4), cell
      real(dp) :: bpos_01(3), bpos_02(3), bpos_23(3), bpos_13(3)

      if (present(position_id)) then
        if (position_id < 1 .or. position_id+2 > VSIZE_RPAR) error stop &
          'mesh_append_rectilinear_mesh - position_id out of bounds'
      end if

      block
        integer :: ntiles(2)
        real(dp) :: mesh_size(2), p3(3)
        real(dp), parameter :: beta=1.5_dp
!TODO verify
! - mesh_size > 0
! - dot_produxt(p1-p0,p2-p0) > 0 (vectors are not parallel)
! - cell_size > 0
! - rel_shift between 0 and 0.5
        mesh_size(1) = sqrt(dot_product(p1-p0, p1-p0))
        mesh_size(2) = sqrt(dot_product(p2-p0, p2-p0))
        base1 = (p1-p0)/mesh_size(1)
        base2 = (p2-p0)/mesh_size(2)
        ntiles(1) = max(1, ceiling(mesh_size(1)/cell_size))
        ntiles(2) = max(1, ceiling(mesh_size(2)/cell_size*(0.75)))
        tile_size = mesh_size/ntiles
        allocate(p_corners(ntiles(1)+1, ntiles(2)+1))
        allocate(p_mids(ntiles(1), ntiles(2)))

        ! points in the corners of rectangulars
        do i=1, size(p_corners,1)
          do j=1, size(p_corners,2)
            p_corners(i,j) = this%add_point( p0 + &
              real(i-1)*tile_size(1)*base1 + real(j-1)*tile_size(2)*base2)
          end do
        end do

        ! "positions" of boundary cells (ghost cells), just for visualization
        p3 = p0 + (p1-p0) + (p2-p0)
        bpos_01 = beta/2.0_dp*(p0 + p1) + (1.0_dp-beta)*p2
        bpos_02 = beta/2.0_dp*(p0 + p2) + (1.0_dp-beta)*p1
        bpos_23 = beta/2.0_dp*(p2 + p3) + (1.0_dp-beta)*p0
        bpos_13 = beta/2.0_dp*(p1 + p3) + (1.0_dp-beta)*p0
      end block

      ! points in the middle of rectangulars
      block
        real(dp) :: shift(2)

        do i=1, size(p_mids,1)
          do j=1, size(p_mids,2)
            associate( &
              c1=>this%points(this%index_from_handle( &
                  p_corners(i,j)))%position, &
              c2=>this%points(this%index_from_handle( &
                  p_corners(i+1,j)))%position, &
              c3=>this%points(this%index_from_handle( &
                  p_corners(i,j+1)))%position, &
              c4=>this%points(this%index_from_handle( &
                  p_corners(i+1,j+1)))%position )
              call random_number(shift)
              shift = (2.0*shift - 1.0) * tile_size * rel_shift

              p_mids(i,j) = this%add_point((c1+c2+c3+c4)/4.0_dp + &
                  base1*shift(1) + base2*shift(2))
            end associate
          end do
        end do
      end block

      ! queues for handles to boundary cells
      call bitems_01%initialize(chunksize=size(transfer(cell,INTEGER_MOLD)))
      call bitems_02%initialize(chunksize=size(transfer(cell,INTEGER_MOLD)))
      call bitems_13%initialize(chunksize=size(transfer(cell,INTEGER_MOLD)))
      call bitems_23%initialize(chunksize=size(transfer(cell,INTEGER_MOLD)))

      ! connect points
      block
        p(4) = p_mids(1,1) ! just some valid point that will be ignored
        do j=1, size(p_mids, 2)
          ! cell to the left of the first middle point
          p(1:3) = [p_mids(1,j), p_corners(1,j), p_corners(1,j+1)]
          cell = this%add_cell(p, position_id)
          call add_boundary_cell(bitems_02, bpos_02, cell)

          ! cells between middle points
          do i=1, size(p_mids,1)-1
            p(1:3) = [p_mids(i,j), p_corners(i+1,j), p_mids(i+1,j)]
            cell = this%add_cell(p, position_id)
            p(1:3) = [p_mids(i,j), p_corners(i+1,j+1), p_mids(i+1,j)]
            cell = this%add_cell(p, position_id)
          end do

          ! cell to the right of the last middle point
          p(1:3) = [p_mids(size(p_mids,1),j), p_corners(size(p_mids,1)+1,j), &
              p_corners(size(p_mids,1)+1,j+1)]
          cell = this%add_cell(p, position_id)
          call add_boundary_cell(bitems_13, bpos_13, cell)

          do i=1, size(p_mids,1)
            ! cell bellow middle point
            p(1:3) = [p_mids(i,j), p_corners(i,j), p_corners(i+1,j)]
            cell = this%add_cell(p, position_id)
            if (j==1) call add_boundary_cell(bitems_01, bpos_01, cell)

            ! cell above middle point
            p(1:3) = [p_mids(i,j), p_corners(i,j+1), p_corners(i+1,j+1)]
            cell = this%add_cell(p, position_id)
            if (j==size(p_mids,2)) call add_boundary_cell(bitems_23, bpos_23, cell)
          end do
        end do
      end block

      ! empty queues and store ghost vertice handles to boundary
      allocate(boundary(bitems_01%size() + bitems_02%size() + bitems_13%size() &
          + bitems_23%size()) )
      call consume_handles(bitems_01, i, 1)
      call consume_handles(bitems_02, i, 2)
      call consume_handles(bitems_23, i, 3)
      call consume_handles(bitems_13, i, 4)
      if (i /= size(boundary)) error stop &
        'mesh_append_rectilinear_mesh - boundary size invalid (internal error'

    contains

      subroutine add_boundary_cell(queue, position, cell0)
        type(queue_t), intent(inout) :: queue
        real(dp), intent(in) :: position(3)
        type(graph_handle_t), intent(in) :: cell0

        integer :: vipar(VSIZE_IPAR), eipar(ESIZE_IPAR)
        real(dp) :: vrpar(VSIZE_RPAR), erpar(ESIZE_RPAR)
        type(graph_handle_t) :: boundary_vertex, edge

        vipar = 0
        eipar = 0
        vrpar = 0.0_dp
        erpar = 0.0_dp
        if (present(position_id)) &
            vrpar(position_id:position_id+2) = position
        boundary_vertex = this%add_vertex(vipar, vrpar)
        edge = this%add_edge( &
          this%cells(this%index_from_handle(cell0))%dual_vertex, &
          boundary_vertex, eipar, erpar)
        call queue%enqueue(transfer(boundary_vertex, INTEGER_MOLD))
      end subroutine add_boundary_cell

      subroutine consume_handles(queue, ind, iside)
        type(queue_t), intent(inout) :: queue
        integer, intent(inout) :: ind
        integer, intent(in) :: iside

        if (iside==1) then
          boundary_offset(iside-1) = 0
          ind = 0
        else
          if (boundary_offset(iside-1) /= ind) error stop &
            'mesh_append_rectilinear_mesh - internal error in consume_handles'
        end if
        do while(.not. queue%empty())
          ind = ind + 1
          boundary(ind) = transfer(queue%dequeue(),boundary(ind))
        end do
        boundary_offset(iside) = ind
      end subroutine consume_handles

    end subroutine mesh_append_rectilinear_mesh


    ! ------------------------------------------------
    ! Transient diffusion/conduction TODO here for now
    ! ------------------------------------------------
    subroutine integrate_pde(this, t_start, t_end, dt_comp, dt_out, &
        position_conductance, position_capacitance, bc_label, x_init, x_out, &
        vmask, emask, vselector, eselector, rtol_l2, rtol_linf)
      class(mesh_t), intent(in) :: this
      real(dp), intent(in) :: t_start, dt_comp, dt_out
      real(dp), intent(inout) :: t_end
      integer, intent(in) :: position_conductance, position_capacitance
      integer, intent(in) :: bc_label(:)
      real(dp), intent(in) :: x_init(:)
      real(dp), intent(out), allocatable :: x_out(:,:)
      logical, intent(in), optional :: vmask(:), emask(:)
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      real(dp), intent(in), optional :: rtol_l2, rtol_linf
!
! TODO Documentation block
!
      integer, parameter :: CG_OK=0, CG_MAXITER=1 ! TODO import from graph_smod_flow
      integer, parameter :: BC_NONE = 0
      logical, allocatable :: vmask0(:), emask0(:), is_external(:)
      integer :: iout, iflag, icomp
      real(dp) :: t, dt, finterpol
      real(dp), allocatable :: x_old(:), x_new(:), tmp(:), diag(:)
      logical :: is_out

      if (this%is_directed()) error stop &
        'integrate_pde - undirected graph required'
      if (position_conductance < 1 .or. position_conductance > ESIZE_RPAR) &
        error stop 'integrate_pde - position_conductance out of bounds'
      if (position_capacitance < 1 .or. position_capacitance > VSIZE_RPAR) &
        error stop 'integrate_pde - position_capacitance out of bounds'
      if (size(bc_label) /= this%nvertices) error stop &
        'integrate_pde - size of bc_label is invalid'
      if (size(x_init) /= this%nvertices) error stop &
        'integrate_pde - size of x_init is invalid'

      call this%build_selection_masks(vmask0, emask0, vmask_provided=vmask, &
          emask_provided=emask, vselector=vselector, eselector=eselector)

      associate (nout => ceiling((t_end-t_start)/dt_out) + 1)
        ! nout >= 2 if t_end > t_start
        allocate(x_out(this%nvertices, nout))
        allocate(x_old(this%nvertices), x_new(this%nvertices))
      end associate
      iout = 1
      icomp = 1
      x_out(:,iout) = x_init
      x_old = x_init
      t = t_start

      allocate(is_external(this%nvertices), source=.true.)
      where (vmask0 .and. bc_label==BC_NONE)
        is_external = .false.
      end where

      block
        integer :: icell, jvertex
        type(cell_geometry_t) :: geom
        allocate(diag(this%nvertices), source=0.0_dp)

        do icell=1, this%ncells
          jvertex = this%index_from_handle(this%cells(icell)%dual_vertex)
          if (jvertex==MAP_NULL) error stop 'integrate_pde - could not find dual vertex'
          geom = this%cells(icell)%geometry(this)
         !diag(jvertex) = this%vertices(jvertex)%rpar(position_capacitance) * &
         !    this%cells(icell)%volume(this) / dt_comp
          diag(jvertex) = this%vertices(jvertex)%rpar(position_capacitance) * &
              geom%volume / dt_comp
        end do
      end block
print *, 'DIAG ',diag

      do
        ! update x
        x_new = x_old
 print *, size(diag), size(x_old), size(is_external)
        call conjugate_gradient(this%graph_t, x_new, position_conductance, &
            is_external, emask0, iflag, diag=diag, x_old=x_old, &
            rtol_l2=rtol_l2, rtol_linf=rtol_linf)
        if (iflag/=CG_OK .and. iflag/=CG_MAXITER) then
          print *, 'conjugate_gradient iflag = ',iflag, t
          error stop 'integration_pde - could not solve step'
        else if (iflag==CG_MAXITER) then
          print *, 'conjugate_gradient WARNING tolerance not met ', t
        end if
        ! update time
        t = t_start + real(icomp,dp)*dt_comp
        icomp = icomp + 1
        ! write output and update "iout"
        if (t >= t_start + real(iout,dp)*dt_out) then
          finterpol = (t - (t_start+real(iout,dp)*dt_out)) / dt_comp
          iout = iout + 1
          x_out(:,iout) = (1.0_dp-finterpol)*x_new + finterpol*x_old
          is_out = .true.
        else
          is_out = .false.
        end if
        ! if end of loop, update "t_end" and write last "x_out" column
        if (t >= t_end) then
          t_end = t
          if (.not. is_out) then
            iout = iout + 1
            if (iout /= size(x_out,2)) then
              print *, iout, shape(x_out)
              error stop 'integrate_pde - internal errro'
            end if
            x_out(:,iout) = x_new
          end if
          exit
        end if
        ! swap x_old with x_new
        call move_alloc(x_old, tmp)
        call move_alloc(x_new, x_old)
        call move_alloc(tmp, x_new)
      end do

      if (iout /= size(x_out,2)) error stop &
        'integrate_pde - not all positions written'

    contains
    end subroutine integrate_pde


    ! ======================
    ! Geometrical primitives
    ! ======================

    pure function vec_product(u, v) result(w)
      real(dp), intent(in) :: u(3), v(3)
      real(dp) :: w(3)
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
    end function vec_product


    pure subroutine solve_3x3(A, b, x, iflag)
      real(dp), intent(in) :: A(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: iflag  ! return zero if ok
!
! Solve a linear system A*x = b for problem size = 3
! Return IFLAG_SINGULAR if pivot is very small indicating a degenerate system.
!
      integer, parameter :: IFLAG_OK = 0, IFLAG_SINGULAR = 1
      real(dp) :: Awrk(3,3), bwrk(3), tmpA(3), tmpb
      integer :: i, j, k, n

      iflag = IFLAG_SINGULAR ! to be able to just return if problem detected
      x = 0.0_dp
      n = size(A,1)
      if (n/=3) error stop 'solve 3x3 - problem size 3 expected'
      if (size(A,2)/=n .or. size(b)/=n .or. size(x)/=n) &
          error stop 'solve_3x3 - array size inconsistent'
      Awrk = A
      bwrk = b
      if (maxval(abs(A)) == 0.0_dp) return

      do i = 1, n
        ! Find pivot position
        k = maxloc(abs(Awrk(i:n,i)),dim=1)+i-1
        ! Swap pivot with i-th row
        if (k /= i) then
          tmpA = Awrk(k,:)
          Awrk(k,:) = Awrk(i,:)
          Awrk(i,:) = tmpA
          tmpb = bwrk(k)
          bwrk(k) = bwrk(i)
          bwrk(i) = tmpb
        end if
        ! Verify regularity
        if (abs(Awrk(i,i)) <= NUMERICAL_TOL) return
        ! Normalize pivot (update b before destroying pivot multiplier)
        bwrk(i) = bwrk(i) / Awrk(i,i)
        Awrk(i,i:n) = Awrk(i,i:n) / Awrk(i,i)
        ! Gaussian elimination of items in i-th column
        do j = i+1, n
          ! bwrk must be updated first
          bwrk(j) = bwrk(j) - Awrk(j,i) * bwrk(i)
          Awrk(j,i:n) = Awrk(j,i:n) - Awrk(j,i) * Awrk(i,i:n)
        end do
      end do
      ! Back substitution
      do i = n, 1, -1
        x(i) = bwrk(i)
        do j = i-1, 1, -1
          bwrk(j) = bwrk(j) - x(i)*Awrk(j,i)
        end do
      end do
      iflag = IFLAG_OK
    end subroutine solve_3x3


    pure function triangle_area_vector(p1, p2, p3) result(avec)
      real(dp), intent(in) :: p1(3), p2(3), p3(3)
      real(dp) :: avec(3)
!
! TODO Documentation
!
      real(dp) :: u(3), v(3)

      u = p2 - p1
      v = p3 - p1
      avec = 0.5_dp * vec_product(u, v)
    end function triangle_area_vector


    pure function circumcentre_circle(p1, p2, p3) result (c)
      real(dp), intent(in) :: p1(3), p2(3), p3(3)
      real(dp) :: c(3)
!
! Centre of a circle with points P1, P2 and P3 on it.
! Points must not be colinear
!
      real(dp) :: u(3), v(3), w(3), unorm, vnorm, wnorm, nom(3)

      u = p2 - p1
      v = p3 - p1
      w = vec_product(u, v)
      unorm = dot_product(u, u)
      vnorm = dot_product(v, v)
      wnorm = dot_product(w, w)
      nom = vnorm*vec_product(w,u) + unorm*vec_product(v,w)
      c = p1 + nom / (2.0 * wnorm)
    end function circumcentre_circle


    pure function circumcentre_sphere(p1, p2, p3, p4) result(c)
      real(dp), intent(in) :: p1(3), p2(3), p3(3), p4(3)
      real(dp) :: c(3)
!
! Center of the sphere with p1, p2, p3 and p4 on its surface.
! Points must not be co-planar.
!
      real(dp) :: u(3), v(3), w(3), A(3,3), b(3), x(3)
      integer :: iflag

      u = p2 - p1
      v = p3 - p1
      w = p4 - p1
      A(1,:) = u
      A(2,:) = v
      A(3,:) = w
      b(1) = 0.5_dp * dot_product(u, u)
      b(2) = 0.5_dp * dot_product(v, v)
      b(3) = 0.5_dp * dot_product(w, w)
      call solve_3x3(A, b, x, iflag)
      if (iflag /= 0) error stop 'circumcentre_sphere - degenerate positions'
      c = p1 + x
    end function circumcentre_sphere


    pure function point_line_distance(p1, p2, p3) result(d)
      real(dp), intent(in) :: p1(3), p2(3), p3(3)
      real(dp) :: d
!
! Shortest distance "d" of P3 from an infinite line through P1 and P2
! P1 and P2 must not be co-incident.
!
! https://paulbourke.net/geometry/pointlineplane/
!
! 1. Assume point P on line P1-P2 that is closest to P3
! 2. P1-P2 and P-P3 are perpendicular, i.e. (P3-P)*(P2-P1) = 0
! 3. P = P1 + u*(P2-P1)
! 4. Substitute [3] to [2]:
!      (P3-P) * (P2-P1) = 0
!      [P3 - P1 - u*(P2-P1)] * (P2-P1) = 0
!      (P3-P1) * (P2-P1) - u*(P2-P1) * (P2-P1) = 0
!          (P3-P1) * (P2-P1)
!      u = -----------------
!             |(P2-P1)|**2
!
      real(dp) :: u, n2, p(3), v(3)

      v = p2-p1
      n2 = dot_product(v, v)
#ifdef DEBUG
      if (n2 <= NUMERICAL_TOL) error stop &
          'point_line_distance - P1 and P2 too close (NUMERICAL_TOL)'
#endif
      u = dot_product(p3-p1, v) / n2
      p = p1 + u*v
      d = norm2(p3-p)
    end function point_line_distance


    pure function point_plane_distance(p1, p2, p3, p4) result(d)
      real(dp), intent(in) :: p1(3), p2(3), p3(3), p4(3)
      real(dp) :: d
!
! Distance "d" of P4 from plane defined by points P1,P2,P3.
! Sign of "d" determines on which side of the plane P4 lies.
! Following P1->P2->P3 with the right-hand fingers, the thumb points
! in the positive direction ("above" the plane).
! Points P1, P2 and P3 must not be co-linear.
!
      real(dp) :: u(3), v(3), n(3), nnorm

      u = p2 - p1
      v = p3 - p1

      n = vec_product(u, v)
      nnorm = dot_product(n, n)
#ifdef DEBUG
      if (nnorm <= NUMERICAL_TOL) &
          error stop 'point_plane_distance - P1,P2,P3 degenerate (NUMERICAL_TOL)'
#endif

      ! Vector "n" points in the normal direction of P1-P2-P3 plane.
      ! Distance of P4 from the plane is the projection of P4-P1 onto
      ! the normalized vector "n"
      d = dot_product(p4-p1, n) / sqrt(nnorm)
    end function point_plane_distance


    pure function points_coincident(p1, p2, tol) result(is_coincident)
      real(dp), intent(in) :: p1(3), p2(3), tol
      logical :: is_coincident
!
! True if P1 and P2 are closer than "tol".
!
      is_coincident = dot_product(p2-p1, p2-p1) <= tol**2
    end function points_coincident


    pure function points_colinear(p1, p2, p3, tol) result(is_colinear)
      real(dp), intent(in) :: p1(3), p2(3), p3(3), tol
      logical :: is_colinear
!
! True if P3 closer to line P1-P2 than "tol".
!
      if (points_coincident(p1,p2,tol)) then
        is_colinear = .true.
      else
        is_colinear = point_line_distance(p1, p2, p3) <= tol
      end if
    end function points_colinear


    pure function points_coplanar(p1, p2, p3, p4, tol) result(is_coplanar)
      real(dp), intent(in) :: p1(3), p2(3), p3(3), p4(3), tol
      logical :: is_coplanar
!
! True if P4 closer than "tol" to P1-P2-P3 plane.
!
      if (points_colinear(p1, p2, p3, tol)) then
        is_coplanar = .true.
      else
        is_coplanar = abs(point_plane_distance(p1, p2, p3, p4)) <= tol
      end if
    end function points_coplanar

  end module mesh_mod
