  module graph_mod
    use iso_fortran_env, only : DP => real64, I1B => int8, I8B => int64, output_unit
    use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR
    use conts_mod, only : queue_t, stack_t, pqueue_t, pqueue_handle_t=>handle_t, &
      PQUEUE_MIN, PQUEUE_MAX
    use graph_adjlist_mod, only : adjlist_t, iterator_t
    implicit none (type, external)
    private

    ! Parametrized derived type (PDT) not working reliably with compilers.
    ! To avoid PDT, array sizes required for the actual implementation
    ! are hardcoded in "graph_user.f90" and imported as ?SIZE_?PAR named
    ! constants
    !
    ! Alternativelly, array sizes can be hardcoded here
    ! integer, parameter :: VSIZE_IPAR=?, VSIZE_RPAR=?, ESIZE_IPAR=?, ESIZE_RPAR=?

    ! Named local constants
    integer, parameter :: DEFAULT_ECAPACITY = 10, DEFAULT_VCAPACITY = 5
    integer, parameter :: MAP_NULL = -1, NOT_INITIALIZED = -1
    integer, parameter :: INTEGER_MOLD(0) = [integer ::]

    integer(I1B), parameter :: VERTEX_HANDLE_TYPE = 1_I1B, &
        EDGE_HANDLE_TYPE = 2_I1B, INVALID_HANDLE_TYPE= 0_I1B

    type, public :: handle_t
      private
      integer :: index_to_map = MAP_NULL
      integer :: version = 1
      integer(i1b) :: handle_type = INVALID_HANDLE_TYPE
    contains
      procedure, private :: handle_eq
      generic :: operator(==) => handle_eq
      procedure :: get_index_to_map => handle_get_index_to_map
    end type handle_t

    type, public :: vertex_t
      integer  :: ipar(VSIZE_IPAR)
      real(dp) :: rpar(VSIZE_RPAR)
      type(adjlist_t) :: ngbs ! list of outgoing edge ids
      type(handle_t) :: handle
    end type vertex_t

    type, public :: edge_t
      type(handle_t) :: src_handle, dst_handle
      integer  :: ipar(ESIZE_IPAR)
      real(dp) :: rpar(ESIZE_RPAR)
      type(handle_t) :: handle
    contains
      procedure :: vertex_indices => edge_vertex_indices
    end type

    type, public :: graph_t
      type(vertex_t), allocatable :: vertices(:)
      type(edge_t), allocatable :: edges(:)
      integer, allocatable, private :: vmap(:), emap(:)
          ! storing position of vertices/edges in "vertices"/"edges" arrays
      integer :: nvertices=NOT_INITIALIZED, nedges
      logical, private :: is_directed_graph=.false.
        ! .true. = directed graph (one-way edges)
        ! .false. = undirected graph (edge direction does not matter)
      type(queue_t), private :: free_vhandles, free_ehandles
    contains
      procedure :: initialize => graph_initialize
      procedure :: add_vertex => graph_add_vertex
      procedure :: add_edge   => graph_add_edge
      procedure :: remove_vertex => graph_remove_vertex
      procedure :: remove_edge => graph_remove_edge
      procedure :: remove_orphaned_edges => graph_remove_orphaned_edges
      procedure :: copy => graph_copy
      procedure :: find_edge_id => graph_find_edge_id
      procedure :: print => graph_print
      procedure :: connected_components => graph_connected_components
      procedure :: shortest_path => graph_shortest_path
      procedure :: maxflow => graph_maxflow
      procedure :: betweenness => graph_betweenness
      procedure :: mincut => graph_mincut
      procedure :: build_selection_masks => graph_build_selection_masks
      procedure :: select_vertices => graph_select_vertices
      procedure :: select_edges => graph_select_edges
      procedure :: is_directed => graph_is_directed
    end type graph_t


    ! vertex and edge selector functions interface
    abstract interface
      pure function is_edge_selected(this, edge) result(is)
        import graph_t, edge_t
        implicit none
        class(graph_t), intent(in) :: this
        type(edge_t), intent(in) :: edge
        logical :: is
      end function

      pure function is_vertex_selected(this, vertex) result(is)
        import graph_t, vertex_t
        implicit none
        class(graph_t), intent(in) :: this
        type(vertex_t), intent(in) :: vertex
        logical :: is
      end function
    end interface

  contains

    ! ---------------------
    ! Handle implementation
    ! ---------------------
    pure function handle_eq(a, b) result(eq)
      class(handle_t), intent(in) :: a, b
      logical :: eq
      eq = a%version==b%version .and. &
          a%index_to_map==b%index_to_map .and. &
          a%handle_type==b%handle_type
    end function handle_eq


    subroutine borrow_handle(this, handle_type, handle)
      class(graph_t), intent(inout) :: this
      integer(i1b), intent(in) :: handle_type
      type(handle_t), intent(out) :: handle

      select case(handle_type)
      case(VERTEX_HANDLE_TYPE)
        if (this%free_vhandles%size()==0) call increase_vertices_capacity(this)
        if (this%free_vhandles%size()==0) error stop 'borrow_handle - no more V-handles available'
        handle = transfer(this%free_vhandles%dequeue(), handle)
      case(EDGE_HANDLE_TYPE)
        if (this%free_ehandles%size()==0) call increase_edges_capacity(this)
        if (this%free_ehandles%size()==0) error stop 'borrow_handle - no more E-handles available'
        handle = transfer(this%free_ehandles%dequeue(), handle)
      case default
        error stop 'borrow_handle: invalid handle_type'
      end select
    end subroutine borrow_handle


    subroutine return_handle(this, handle)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: handle

      type(handle_t) :: reused_handle

      reused_handle = handle
      reused_handle%version = reused_handle%version + 1
      select case(reused_handle%handle_type)
      case(VERTEX_HANDLE_TYPE)
        call this%free_vhandles%enqueue(transfer(reused_handle,INTEGER_MOLD))
      case(EDGE_HANDLE_TYPE)
        call this%free_ehandles%enqueue(transfer(reused_handle,INTEGER_MOLD))
      case default
        error stop 'return_handle: handle provided has invalid handle_type'
      end select
    end subroutine return_handle


    pure function get_index_from_handle(this, handle) result(id)
      class(graph_t), intent(in) :: this
      type(handle_t), intent(in) :: handle
      integer id
!
! Return position of a vertex/edge in array using handle. If handle refers to
! the vertex/edge that is no longer in array, MAP_NULL is returned.
!
      id = MAP_NULL
      select case(handle%handle_type)
      case(VERTEX_HANDLE_TYPE)
        if (handle%index_to_map > 0 .and. handle%index_to_map <= size(this%vmap)) then
          id = this%vmap(handle%index_to_map)
          if (id/=MAP_NULL) then
            ! verify version matches the stored one
            if (.not. (this%vertices(id)%handle==handle)) id = MAP_NULL
          end if
        end if
      case(EDGE_HANDLE_TYPE)
        if (handle%index_to_map > 0 .and. handle%index_to_map <= size(this%emap)) then
          id = this%emap(handle%index_to_map)
          if (id/=MAP_NULL) then
            ! verify version matches the stored one
            if (.not. (this%edges(id)%handle==handle)) id = MAP_NULL
          end if
        end if
      case default
        error stop 'get_index_from_handle: unknown handle_type'
      end select
    end function get_index_from_handle


    elemental integer function handle_get_index_to_map(this) result(id)
      class(handle_t), intent(in) :: this
      id = this%index_to_map
    end function


    ! ----------------------
    ! Graph basic operations
    ! ----------------------
    subroutine graph_initialize(this, vcapacity, ecapacity, is_directed_graph)
      class(graph_t), intent(inout) :: this
      integer, intent(in), optional :: vcapacity, ecapacity
      logical, intent(in), optional :: is_directed_graph

      ! directed or undirected gtaph?
      this%is_directed_graph = .false.
      if (present(is_directed_graph)) this%is_directed_graph = is_directed_graph

      ! reallocate all arrays to zero size
      if (allocated(this%vertices)) deallocate(this%vertices)
      allocate(this%vertices(0))
      if (allocated(this%edges)) deallocate(this%edges)
      allocate(this%edges(0))
      if (allocated(this%vmap)) deallocate(this%vmap)
      allocate(this%vmap(0))
      if (allocated(this%emap)) deallocate(this%emap)
      allocate(this%emap(0))

      this%nvertices = 0
      this%nedges = 0

      block ! initialize queues of free handles
        type(handle_t) :: handle
        call this%free_vhandles%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
        call this%free_ehandles%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
      end block

      block ! allocate initial capacity
        integer :: new_capacity
        new_capacity = DEFAULT_VCAPACITY
        if (present(vcapacity)) new_capacity = vcapacity
        call increase_vertices_capacity(this, new_capacity)
        new_capacity = DEFAULT_ECAPACITY
        if (present(ecapacity)) new_capacity = ecapacity
        call increase_edges_capacity(this, new_capacity)
      end block
    end subroutine graph_initialize


    pure logical function graph_is_initialized(this)
      class(graph_t), intent(in) :: this
      graph_is_initialized = this%nvertices /= NOT_INITIALIZED
    end function graph_is_initialized


    pure logical function graph_is_directed(this) result(is)
      class(graph_t), intent(in) :: this
!
! Return TRUE if directerd graph, return FALSE if undirected graph.
! This is a getter function - to be used outside module
!
      if (graph_is_initialized(this)) then
        is = this%is_directed_graph
      else
        error stop 'graph_is_directed - graph not initialized'
      end if
    end function graph_is_directed


    subroutine increase_vertices_capacity(this, new_capacity)
      class(graph_t), intent(inout) :: this
      integer, intent(in), optional :: new_capacity

      integer :: old_capacity, new_capacity0
      type(vertex_t), allocatable :: tmp_vertices(:)
      integer, allocatable :: tmp_map(:)

      old_capacity = size(this%vertices)
      if (present(new_capacity)) then
        new_capacity0 = new_capacity
      else
        new_capacity0 = 2*old_capacity
      end if

      ! reallocate "vertices" and "vmap"
      allocate(tmp_vertices(new_capacity0))
      allocate(tmp_map(new_capacity0), source=MAP_NULL)
      tmp_vertices(1:old_capacity) = this%vertices
      tmp_map(1:old_capacity) = this%vmap
      call move_alloc(tmp_vertices, this%vertices)
      call move_alloc(tmp_map, this%vmap)

      block ! create fresh handles
        integer :: i
        type(handle_t) :: new_handle
        do i=old_capacity+1, new_capacity0
          new_handle%index_to_map = i
          new_handle%version = 1
          new_handle%handle_type = VERTEX_HANDLE_TYPE
          call this%free_vhandles%enqueue(transfer(new_handle, INTEGER_MOLD))
        end do
      end block
    end subroutine increase_vertices_capacity


    subroutine increase_edges_capacity(this, new_capacity)
      class(graph_t), intent(inout) :: this
      integer, intent(in), optional :: new_capacity

      integer :: old_capacity, new_capacity0
      type(edge_t), allocatable :: tmp_edges(:)
      integer, allocatable :: tmp_map(:)

      old_capacity = size(this%edges)
      if (present(new_capacity)) then
        new_capacity0 = new_capacity
      else
        new_capacity0 = 2*old_capacity
      end if

      ! reallocate "edges" and "vmap"
      allocate(tmp_edges(new_capacity0))
      allocate(tmp_map(new_capacity0), source=MAP_NULL)
      tmp_edges(1:old_capacity) = this%edges
      tmp_map(1:old_capacity) = this%emap
      call move_alloc(tmp_edges, this%edges)
      call move_alloc(tmp_map, this%emap)

      block ! create fresh handles
        integer :: i
        type(handle_t) :: new_handle
        do i=old_capacity+1, new_capacity0
          new_handle%index_to_map = i
          new_handle%version = 1
          new_handle%handle_type = EDGE_HANDLE_TYPE
          call this%free_ehandles%enqueue(transfer(new_handle, INTEGER_MOLD))
        end do
      end block
    end subroutine increase_edges_capacity


    function graph_add_vertex(this, ipar, rpar) result(handle)
      class(graph_t), intent(inout) :: this
      integer, intent(in) :: ipar(:)
      real(dp), intent(in) :: rpar(:)
      type(handle_t) :: handle

      if (.not. graph_is_initialized(this)) then
        error stop 'graph_add_vertex - graph not initialized'
      else if (size(ipar) /= VSIZE_IPAR .or. size(rpar) /= VSIZE_RPAR) then
        error stop 'graph_add_vertex - invalid argument arrays size'
      end if

      call borrow_handle(this, VERTEX_HANDLE_TYPE, handle)
      this%nvertices = this%nvertices + 1
      associate (new_vertex => this%vertices(this%nvertices))
        new_vertex%ipar = ipar
        new_vertex%rpar = rpar
        call new_vertex%ngbs%initialize()
        new_vertex%handle = handle
      end associate
      this%vmap(handle%index_to_map) = this%nvertices
    end function graph_add_vertex


    subroutine graph_remove_vertex(this, handle)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: handle
!
! Remove vertex. Edges associated with the vertex will be also removed.
! In directed graphs, only outgoing edges will be removed and the incomming
! edges become orphaned and must be removed by a separate operation.
!
      integer :: ivertex

      if (handle%handle_type /= VERTEX_HANDLE_TYPE) &
          error stop 'graph_remove_vertex - invalid handle type'
      ivertex = get_index_from_handle(this, handle)
      if (ivertex == MAP_NULL) &
          error stop 'graph_remove_vertex - vertex no longer present in graph'

      ! Automatically remove all outgoing edges.
      ! For directed graphs, the incomming edges will become orphaned and must
      ! be removed manually.
      block
        type(iterator_t) :: iterator
        integer :: iedge
        do
          iterator = iterator_t()
          if (.not. this%vertices(ivertex)%ngbs%has_next(iterator)) exit
          call this%vertices(ivertex)%ngbs%next(iterator, iedge)
          call graph_remove_edge(this, this%edges(iedge)%handle)
        end do
      end block

      ! Defensive - could be removed later
      if (this%vertices(ivertex)%ngbs%size()>0) &
          error stop 'graph_remove_vertex - could not remove outgoing edges'

      ! Nullify vmap and return handle
      this%vmap(handle%index_to_map) = MAP_NULL
      call return_handle(this, handle)

      ! Relocate the last vertex to fill "hole" after removed vertex
      if (ivertex /= this%nvertices) then
        call relocate_vertex(this, this%vertices(this%nvertices)%handle, ivertex)
      end if
      this%nvertices = this%nvertices - 1
    end subroutine graph_remove_vertex


    subroutine relocate_vertex(this, handle, newid)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: handle
      integer, intent(in) :: newid

      integer :: oldid

      if (handle%handle_type /= VERTEX_HANDLE_TYPE) &
          error stop 'relocate_vertex - wrong handle type'
      oldid = get_index_from_handle(this, handle)
      if (oldid == MAP_NULL) &
          error stop 'relocate_vertex - vertex no more exists'
      if (newid < 1 .or. newid > this%nvertices) &
          error stop 'relocate_vertex - newid out of bounds'

      ! copy vertex and update record in "vmap"
      this%vertices(newid) = this%vertices(oldid)
      this%vmap(handle%index_to_map) = newid
    end subroutine relocate_vertex


    function graph_add_edge(this, src, dst, ipar, rpar) result(handle)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: src, dst
      integer, intent(in) :: ipar(:)
      real(dp), intent(in) :: rpar(:)
      type(handle_t) :: handle

      if (.not. graph_is_initialized(this)) then
        error stop 'graph_add_edge - graph not initialized'
      else if (size(ipar) /= ESIZE_IPAR .or. size(rpar) /= ESIZE_RPAR) then
        error stop 'graph_add_edge - invalid argument arrays size'
      else if (src%handle_type /= VERTEX_HANDLE_TYPE .or. dst%handle_type /= VERTEX_HANDLE_TYPE) then
        error stop 'graph_add_edge - invalid src or dst handles type'
      end if

      block
        integer :: isrc, idst

        isrc = get_index_from_handle(this, src)
        idst = get_index_from_handle(this, dst)
        if (isrc==MAP_NULL .or. idst==MAP_NULL) then
          error stop 'graph_add_edge - vertex not present (invalid handle)'
        end if
        ! check if connection already exists
        if (this%find_edge_id(isrc, idst) /= MAP_NULL) then
          error stop 'graph_add_edge - src-dst connection already exists'
        else if (.not. this%is_directed_graph) then
          ! undirected graph: look also in destination's vertex ngb-list
          if (this%find_edge_id(idst, isrc) /= MAP_NULL) then
            error stop 'graph_add_edge - opposite connection already exists'
          end if
        end if

        ! we verified there is no connection between SRC and DST
        call borrow_handle(this, EDGE_HANDLE_TYPE, handle)
        this%nedges = this%nedges + 1
        associate (new_edge => this%edges(this%nedges))
          new_edge%src_handle = src
          new_edge%dst_handle = dst
          new_edge%ipar = ipar
          new_edge%rpar = rpar
          new_edge%handle = handle
        end associate
        this%emap(handle%index_to_map) = this%nedges

        ! add the new edge to the source vertex's outgoing edges list...
        call this%vertices(isrc)%ngbs%add(this%nedges)
        if (.not. this%is_directed_graph) then
          ! ...and also to the destination vertex's outgoing edges list
          ! for an undirected graph
          call this%vertices(idst)%ngbs%add(this%nedges)
        end if
      end block
    end function graph_add_edge


    subroutine graph_remove_edge(this, handle)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: handle

      integer :: iedge, isrc, idst

      if (handle%handle_type /= EDGE_HANDLE_TYPE) &
          error stop 'graph_remove_edge - invalid handle type'
      iedge = get_index_from_handle(this, handle)
      if (iedge == MAP_NULL) &
          error stop 'graph_remove_edge - edge no longer exists'

      ! Remove edge refrence from adjacency list(s) of end points
      ! It is ok if the vertice no longer exist, but if it exists, the
      ! reference to the edge must be present.
      isrc = get_index_from_handle(this, this%edges(iedge)%src_handle)
      idst = get_index_from_handle(this, this%edges(iedge)%dst_handle)
      if (isrc /= MAP_NULL) call this%vertices(isrc)%ngbs%remove(iedge)
      if (.not. this%is_directed_graph) then
        if (idst /= MAP_NULL) call this%vertices(idst)%ngbs%remove(iedge)
      end if

      ! Nullify emap entry and reuse the handle
      this%emap(handle%index_to_map) = MAP_NULL
      call return_handle(this, handle)

      ! Relocate the last edge to the "hole" after removed edge
      if (iedge /= this%nedges) then
        call relocate_edge(this, this%edges(this%nedges)%handle, iedge)
      end if
      this%nedges = this%nedges - 1
    end subroutine graph_remove_edge


    subroutine relocate_edge(this, handle, newid)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: handle
      integer, intent(in) :: newid

      integer :: oldid, isrc, idst

      if (handle%handle_type /= EDGE_HANDLE_TYPE) &
          error stop 'relocate_edge - wrong handle type'
      oldid = get_index_from_handle(this, handle)
      if (oldid == MAP_NULL) &
          error stop 'relocate_edge - edge no more exists'
      if (newid < 1 .or. newid > this%nedges) &
          error stop 'relocate_edge - newid out of bounds'

      ! copy edge and update record in "emap"
      this%edges(newid) = this%edges(oldid)
      this%emap(handle%index_to_map) = newid

      ! update adjacent lists of respective vertices
      isrc = get_index_from_handle(this, this%edges(newid)%src_handle)
      if (isrc/=MAP_NULL) call update_ngbs(this%vertices(isrc)%ngbs)

      if (.not. this%is_directed_graph) then
        idst = get_index_from_handle(this, this%edges(newid)%dst_handle)
        if (idst/=MAP_NULL) call update_ngbs(this%vertices(idst)%ngbs)
      end if

    contains
      subroutine update_ngbs(ngbs) ! internal procedure
        type(adjlist_t), intent(inout) :: ngbs

        type(iterator_t) :: found_oldid
        if (ngbs%contains(newid)) &
            error stop 'relocate edge - newid present in list would lead to duplicity'
        found_oldid = ngbs%find(oldid)
        if (.not. ngbs%has_next(found_oldid)) &
            error stop 'relocate edge - old id not found in adjacent list'
        call ngbs%remove(oldid, found_oldid)
        call ngbs%add(newid, skip_duplicity_check=.true.)
      end subroutine
    end subroutine relocate_edge


    subroutine graph_remove_orphaned_edges(this, nedges_removed)
      class(graph_t), intent(inout) :: this
      integer, intent(out), optional :: nedges_removed
!
! This subroutine must be called after removing vertices in a directed graph
! to remove edges whose destination vertex is no longer present.
! Orphaned edges should not be present in undireceted graph.
!
      integer :: iedge, isrc, idst, nedges_removed0

      nedges_removed0 = 0
      do iedge=1, this%nedges
        isrc = get_index_from_handle(this, this%edges(iedge)%src_handle)
        idst = get_index_from_handle(this, this%edges(iedge)%dst_handle)
        if (isrc==MAP_NULL) error stop &
          'remove_orphaned_edges - non-existing source vertex is unexpected'
        if (idst==MAP_NULL) then
          if (.not. this%is_directed_graph) error stop &
            'remove_orphaned_edges - non-existing destination vertex is unexpected in undirected graphs'
          call this%remove_edge(this%edges(iedge)%handle)
          nedges_removed0 = nedges_removed0 + 1
        end if
      end do
      if (present(nedges_removed)) nedges_removed = nedges_removed0
print '("temove_orphaned_edges: removed ",i0," edges")', nedges_removed0
    end subroutine graph_remove_orphaned_edges


    subroutine graph_copy(this, gnew, vselector, eselector, vmask, emask, &
        new_vertices, new_edges)
      class(graph_t), intent(in) :: this
      type(graph_t), intent(inout) :: gnew
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      logical, intent(in), optional :: vmask(:), emask(:)
      type(handle_t), intent(inout), allocatable, optional :: &
          new_vertices(:), new_edges(:)
!
! Copy selected vertices and edges to a new graph
!
      logical, allocatable :: vmask0(:), emask0(:)
      type(handle_t), allocatable :: new_vertices0(:), new_edges0(:)

      ! Select copied vertices and edges
      call this%build_selection_masks(vmask0, emask0, &
          vmask_provided=vmask, emask_provided=emask, &
          vselector=vselector, eselector=eselector)

      ! Initialize an empty graph and arrays maping vertices/edges in the 
      ! new graph to their counterparts in the old graph.
      call gnew%initialize(is_directed_graph=this%is_directed_graph, &
          vcapacity=count(vmask0), ecapacity=count(emask0))
      allocate(new_vertices0(this%nvertices), source=handle_t(MAP_NULL,0,VERTEX_HANDLE_TYPE))
      allocate(new_edges0(this%nedges), source=handle_t(MAP_NULL,0,EDGE_HANDLE_TYPE))

      ! Copy selected vertices
      block
        integer :: v
        do v=1, this%nvertices
          if (.not. vmask0(v)) cycle
          associate(vertex=>this%vertices(v))
            new_vertices0(v)=gnew%add_vertex(vertex%ipar, vertex%rpar)
          end associate
        end do
      end block

      ! Copy selected edges
      block
        integer :: e, ia, ib
        do e=1, this%nedges
          if (.not. emask0(e)) cycle
          associate(edge=>this%edges(e))
            ia = get_index_from_handle(this, edge%src_handle)
            ib = get_index_From_handle(this, edge%dst_handle)
            new_edges0(e) = gnew%add_edge( &
                new_vertices0(ia), new_vertices0(ib), edge%ipar, edge%rpar)
          end associate
        end do
      end block

      ! Return handles to vertices/edges in new graph
      if (present(new_vertices)) call move_alloc(new_vertices0, new_vertices)
      if (present(new_edges)) call move_alloc(new_edges0, new_edges)

    end subroutine graph_copy


    subroutine graph_print(this, fid)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: fid
!
! Print all graph data for debugging
!
      integer :: i, j, v1, v2
      integer, allocatable :: ngbsid(:)
      character(len=:), allocatable :: str_graph_type

      if (this%is_directed_graph) then
        str_graph_type = 'Directed graph'
      else
        str_graph_type = 'Undirected graph'
      end if
      write(fid,'("--- graph dump ---")')
      write(fid, '(a," with ",i0," vertices and ",i0," edges")') &
        str_graph_type, this%nvertices, this%nedges

      ! information about vertices
      do i=1, this%nvertices
        write(fid, '("V-",i0,", connected to")', advance='no') &
          this%vertices(i)%handle%index_to_map
        if (allocated(ngbsid)) deallocate(ngbsid)
        allocate(ngbsid(this%vertices(i)%ngbs%size()))
        ngbsid = list_of_ngbs(this, i)
        do j=1, size(ngbsid)
          write(fid,'(" V-",i0)', advance='no') &
            this%vertices(ngbsid(j))%handle%index_to_map
        end do
        if (size(ngbsid)==0) then
          write(fid,'(" no one:")')
        else
          write(fid,'(":")')
        end if
        if (VSIZE_IPAR>0) write(fid,*) this%vertices(i)%ipar
        if (VSIZE_RPAR>0) write(fid,*) this%vertices(i)%rpar
      end do

      ! information about edges
      do i=1, this%nedges
        v1 = get_index_from_handle(this, this%edges(i)%src_handle)
        v2 = get_index_from_handle(this, this%edges(i)%dst_handle)
        if (v1==MAP_NULL) then
          v1=v1
        else
          v1=this%vertices(v1)%handle%index_to_map
        end if
        if (v2==MAP_NULL) then
          v2=v2
        else
          v2=this%vertices(v2)%handle%index_to_map
        end if
        write(fid, '("E-",i0," connecting V-",i0," and V-",i0,":")') &
          this%edges(i)%handle%index_to_map, v1, v2
        if (ESIZE_IPAR>0) write(fid,*) this%edges(i)%ipar
        if (VSIZE_RPAR>0) write(fid,*) this%edges(i)%rpar
      end do

      write(fid,'("--- end of graph dump ---",/)')
    end subroutine graph_print


    ! ----------------------
    ! Graph helper functions
    ! ----------------------
    pure function list_of_ngbs(this, isrc) result(idsts)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: isrc
      integer :: idsts(this%vertices(isrc)%ngbs%size())
!
! Return an array of neighbors of "isrc" vertex.
!
      integer :: ipos, iedge, idst
      type(iterator_t) :: iterator

      iterator = iterator_t()
      ipos = 0
      do while (this%vertices(isrc)%ngbs%has_next(iterator))
        call this%vertices(isrc)%ngbs%next(iterator, iedge)

        if (iedge <= 0 .or. iedge > this%nedges) then
          error stop 'list_of_ngbs - item in ngbs is out of bounds'
        end if

        idst = other_vertex_id(this, iedge, isrc)
        ! for directed graphs verify, that destination vertex has been selected
        if (this%is_directed_graph) then
          if (idst /= get_index_from_handle(this, this%edges(iedge)%dst_handle)) &
            error stop 'list_of_ngbs - wrong source in directed graph'
        end if

        ipos = ipos+1
        if (ipos > size(idsts)) error stop 'list_of_ngbs - something wrong'
        idsts(ipos) = idst
      end do
    end function list_of_ngbs


    pure function list_of_outgoing_edges(this, isrc) result(iedges)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: isrc
      integer :: iedges(this%vertices(isrc)%ngbs%size())
!
! Return an array of edge positions in "edges" array.
!
      type(iterator_t) :: iterator
      integer :: i

      iterator = iterator_t()
      do i=1, this%vertices(isrc)%ngbs%size()
        if (.not. this%vertices(isrc)%ngbs%has_next(iterator)) &
            error stop 'list_of_outgoing_edges - something wrong'
        call this%vertices(isrc)%ngbs%next(iterator, iedges(i))
      end do
    end function list_of_outgoing_edges


    pure function other_vertex_id(this, iedge, ia, allow_orphaned_edge) result(ib)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: iedge, ia
      integer ib
      logical, intent(in), optional :: allow_orphaned_edge !default = .false.
!
! Given the edge and one of its end-point vertices, return index of the other
! end-point vertex.
! Throws an error if other end-point vertex no longer exists unless optional
! "allow_orphaned_edge" argument is set to .true.
!
! INPUT
!   this  - graph object
!   iedge - position of the edge in "edges" array
!   ia    - position of one vertex in "vertices" array
!   allow_orphaned_edge - (optional) other vertex may no longer exist
! OUTPUT
!   ib    - position of the other vertex in "vertices" array
!
      integer :: i1, i2
      logical :: ignore

      if (ia==MAP_NULL) &
          error stop 'other_vertex_id - "ia" must not be null'

      i1 = get_index_from_handle(this, this%edges(iedge)%src_handle)
      i2 = get_index_from_handle(this, this%edges(iedge)%dst_handle)
      if (ia==i1) then
        ib = i2
      else if (ia==i2) then
        ib = i1
      else
        error stop 'other_vertex_id - no edge endpoint vertex matches to "ia"'
      end if

      if (ib==MAP_NULL) then
        ! How an orphaned edge is treated?
        ignore = .false.
        if (present(allow_orphaned_edge)) ignore = allow_orphaned_edge
        if (.not. ignore) error stop &
            'other_vertex_id - other end point vertex no longer exists'
      end if
    end function other_vertex_id


    pure function graph_find_edge_id(this, ia, ib) result(id)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: ia, ib
      integer :: id
!
! Given two vertices, find the edge that connects them and return its index.
! Return "MAP_NULL" if no such edge exists.
!
! Only an edge from "ia" to "ib" is found in directed graph.
!
! INPUT
!   this   - graph object
!   ia, ib - position of two vertices in "vertices" array
! OUTPUT
!   id     - position of the found edge in "edges" array or MAP_NULL
!
      integer :: iedge
      type(iterator_t) :: iterator

      id = MAP_NULL
      ! Iterate through adjacency list of "ia" vertex
      iterator = iterator_t()
      do while (this%vertices(ia)%ngbs%has_next(iterator))
        call this%vertices(ia)%ngbs%next(iterator, iedge)
        if (other_vertex_id(this, iedge, ia) /= ib) cycle
        ! Other vertex is "ib"
        id = iedge
        exit
      end do
    end function graph_find_edge_id


    ! -------------------
    ! Edge TPB procedures
    ! -------------------
    pure function edge_vertex_indices(this, graph) result(ids)
      class(edge_t), intent(in) :: this
      type(graph_t), intent(in) :: graph
      integer ::ids(2)

      ids(1) = get_index_from_handle(graph, this%src_handle)
      ids(2) = get_index_from_handle(graph, this%dst_handle)
    end function edge_vertex_indices


    ! --------------------------
    ! Label connected components
    ! --------------------------
    subroutine graph_connected_components(this, position_label, &
        vselector, eselector, vmask, emask, lab_count)
      class(graph_t), intent(inout) :: this
      integer, intent(in) :: position_label
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      logical, intent(in), optional :: vmask(:), emask(:)
      integer, intent(out), optional :: lab_count
!
! Identify connected components of an undirected graph.
! (NOTE - directed graphs are left as the future work)
!
! The routine assigns an integer component label to each selected vertex.
! Vertices and edges can be restricted using optional selector functions.
! If no selectors are provided, all vertices and edges are considered.
!
! Only selected vertices are modified. Vertices excluded by the selection
! remain unchanged. An edge contributes to connectivity only if it is selected
! and both of its endpoint vertices are selected.
!
! Component labels are assigned consecutively starting from 1.
! Labels of unselected vertices are not modified.
!
! INPUT
!   this           - graph structure
!   position_label - position in the vertices(:)/ipar array where component
!                    labels are stored
!   vselector      - optional function selecting vertices that participate in
!                    the component search
!   eselector      - optional function selecting edges that participate in
!                    the component search
!   vmask, emask   - optional array selecting vetices and edges participating in
!                    the component search
!
! OUTPUT
!   lab_count      - optional number of connected components identified
!
      integer, parameter :: LAB_INPROGRESS=-1
      integer :: i, j, k, iedge, lab_current
      type(iterator_t) :: iterator
      type(stack_t) :: stack
      logical, allocatable :: vmask0(:), emask0(:)

      ! At the moment only for undirected graphs
      if (this%is_directed_graph) &
          error stop 'graph_connected_components - graph must be undirected'

      call graph_build_selection_masks(this, vmask0, emask0, &
          vselector=vselector, eselector=eselector, vmask_provided=vmask, &
          emask_provided=emask)

      ! Initialize "labels" for selected vertices
      where (vmask0) &
          this%vertices(1:this%nvertices)%ipar(position_label) = LAB_INPROGRESS

      ! Identified components counter
      lab_current = 0

      ! Stack for deep-first graph traversal (DFS)
      call stack%initialize(chunksize=size(transfer(i,INTEGER_MOLD)))

      MAIN_LOOP: do i=1, this%nvertices
        ! Find the next unprocessed vertex and add it to the empty stack
        if (.not. vmask0(i)) cycle
        if (this%vertices(i)%ipar(position_label) /= LAB_INPROGRESS) cycle
        lab_current = lab_current + 1
        this%vertices(i)%ipar(position_label) = lab_current
        call stack%push(transfer(i,INTEGER_MOLD))

        ! Process the stack and propagate "lab_current"
        DFS_LOOP: do while (.not. stack%empty())
          j = transfer(stack%pop(), j)

          ! Label and add allowed neighbours to the stack
          iterator = iterator_t()
          NGB_LOOP: do while (this%vertices(j)%ngbs%has_next(iterator))
            call this%vertices(j)%ngbs%next(iterator, iedge)
            if (.not. emask0(iedge)) cycle
            k = other_vertex_id(this, iedge, j)
            ! Assert "k" is selected as edges to unselected vertices should
            ! have been unselected by "build_selection_masks"
            if (.not. vmask0(k)) error stop &
                'graph_connected_components - selected edge has an unselected endpoint'
            ! According to its label, a neighbor "k" can be:
            !   - unvisited selected vertex: assign current component
            !     and push to stack
            !   - already assigned to this component: ignore
            ! Other labels indicate an inconsistent graph traversal.
            associate (lab_dst=>this%vertices(k)%ipar(position_label))
              if (lab_dst==LAB_INPROGRESS) then
                lab_dst = lab_current
                call stack%push(transfer(k,INTEGER_MOLD))
              else if (lab_dst==lab_current) then
                continue
              else
                ! assertion may fail in directed graphs, but should not happen
                ! in undirected graphs
                error stop 'graph_connected_components - neighbour is already labeled, traversal inconsistency'
              end if
            end associate
          end do NGB_LOOP

        end do DFS_LOOP

      end do MAIN_LOOP

      if (present(lab_count)) lab_count = lab_current

    end subroutine graph_connected_components


    ! -----------------------------
    ! Dijkstra shortest path search
    ! -----------------------------
    subroutine graph_shortest_path(this, position_distance, position_cost, &
        start_vertex, target_vertex, vselector, eselector, path)
      class(graph_t), intent(inout) :: this
      integer, intent(in) :: position_distance, position_cost
      type(handle_t), intent(in) :: start_vertex
      type(handle_t), intent(in), optional :: target_vertex
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      type(handle_t), allocatable, intent(out), optional :: path(:)
!
! Find the shortest path between two vertices using Dijkstra's algorithm.
!
! INPUT
!   this              - graph structure, distance field updated on return
!   position_distance - position in "vertices%rpar" for the distance from
!                       the starting vertex on return (OUT)
!   position_cost     - position in "edges%rpar" where the cost traversing
!                       the edge is stored (IN)
!   start_vertex      - handle to the starting vertex
!   target_vertex     - handle to the target vertex (optional)
!   vselector         - user function to control which vertices can be passed
!                       through (optional)
!   eselector         - user function to control which edges can be used
!                       (optional)
! OUTPUT
!   path              - array of vertex handles on the shortest path (optional)
!
      logical, allocatable :: visited(:), vmask(:), emask(:)
      integer, allocatable :: prev_id(:)
      real(dp) :: cost_to_ngb
      integer :: id_start, id_target, i, id_current, id_ngb, iedge
      type(pqueue_t) :: pqueue
      type(pqueue_handle_t), allocatable :: handles(:)
      type(iterator_t) :: iterator

      ! Find starting and target vertex positions in "vertices" array
      id_start = get_index_from_handle(this, start_vertex)
      if (id_start==MAP_NULL .or. start_vertex%handle_type/=VERTEX_HANDLE_TYPE) &
          error stop 'graph_shortest_path - starting vertex not identified'
      if (present(target_vertex)) then
        id_target = get_index_from_handle(this, target_vertex)
        if (id_target==MAP_NULL .or. target_vertex%handle_type/=VERTEX_HANDLE_TYPE) &
            error stop 'graph_shortest_path - target vertex not identified'
      else
        ! target vertex not provided: search shortest path to every vertex
        ! reachable from start
        id_target = MAP_NULL
      end if

      ! Build selection masks
      call graph_build_selection_masks(this, vmask, emask, vselector=vselector, eselector=eselector)

      ! Local working arrays. Set initial values.
      allocate(visited(this%nvertices), source=.false.)
      allocate(prev_id(this%nvertices), source=MAP_NULL)
      this%vertices(1:this%nvertices)%rpar(position_distance) = huge(cost_to_ngb)
      call pqueue%initialize(chunksize=size(transfer(i,INTEGER_MOLD)), ordering=PQUEUE_MIN)
      allocate(handles(this%nvertices))

      ! Insert starting vertex to the queue
      if (vmask(id_start)) then
        associate(d=>this%vertices(id_start)%rpar(position_distance))
          d = 0.0_dp
          handles(id_start) = pqueue%insert(transfer(id_start,INTEGER_MOLD), d)
        end associate
      end if

      MAIN_LOOP: do
        if (pqueue%empty()) exit MAIN_LOOP
        id_current = transfer(pqueue%pop(), i)
        visited(id_current) = .true.
        if (id_current==id_target) exit MAIN_LOOP

        iterator = iterator_t()
        NGB_LOOP: do while (this%vertices(id_current)%ngbs%has_next(iterator))
          call this%vertices(id_current)%ngbs%next(iterator, iedge)
          if (.not. emask(iedge)) cycle
          id_ngb = other_vertex_id(this, iedge, id_current)
          if (.not. vmask(id_ngb)) cycle
          if (visited(id_ngb)) cycle

          ! id_ngb is an unvisited neighbour of id_current
          associate(d=>this%vertices(id_ngb)%rpar(position_distance))
            cost_to_ngb = this%vertices(id_current)%rpar(position_distance) + &
                        & this%edges(iedge)%rpar(position_cost)
            if (cost_to_ngb < d) then
              ! shorter path found to ngb
              d = cost_to_ngb
              prev_id(id_ngb) = id_current
              if (pqueue%contains(handles(id_ngb))) then
                call pqueue%update_priority(handles(id_ngb), cost_to_ngb)
              else
                handles(id_ngb) = pqueue%insert(transfer(id_ngb,INTEGER_MOLD), cost_to_ngb)
              end if
            end if
          end associate
        end do NGB_LOOP
      end do MAIN_LOOP

      ! Back track from target to construct the path
      if (present(path)) then
        if (.not. present(target_vertex)) &
            error stop 'graph_shortest_path - can not return path if target not given'
        block
          type(stack_t) :: stack

          call stack%initialize(chunksize=size(transfer(i,INTEGER_MOLD)))
          id_current = id_target
          do while (prev_id(id_current)/=MAP_NULL)
            call stack%push(transfer(id_current,INTEGER_MOLD))
            id_current = prev_id(id_current)
          end do
          ! add starting vertex to the path
          if (id_current==id_start) call stack%push(transfer(id_current,INTEGER_MOLD))
          ! the stack will be empty, if target vertex is unreachable

          ! copy stack content to the output array
          allocate(path(stack%size()))
          do i=1, size(path)
            if (stack%empty()) error stop 'graph_shortest_path - stack unexpectedly empty'
            path(i) = this%vertices(transfer(stack%pop(),i))%handle
          end do
        end block
      end if

    end subroutine graph_shortest_path


    ! ---------------------------------------
    ! Edmonds-Karp algorithm for maximum flow
    ! ---------------------------------------
    subroutine graph_maxflow(this, source, sink, position_capacity, flow, &
        position_mincutlabel, position_flow, vselector, eselector, vmask, emask)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: source, sink
      integer, intent(in) :: position_capacity
      real(dp), intent(out) :: flow
      integer, intent(in), optional :: position_mincutlabel
      integer, intent(in), optional :: position_flow
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      logical, intent(in), optional :: vmask(:), emask(:)
!
! Maximum flow from the source to sink.
!
! For directed graphs, temporary reverse edges are added internally to
! represent the residual network. These edges are removed before the
! routine returns.
!
! INPUT
!   this              - the graph (vertex/edge data updated)
!   source            - handle to the source vertex
!   sink              - handle to the sink vertex
!   position_capacity - "edges/rpar" array item giving the edge capacity
!   vselector         - user function to select open verices (OPTIONAL)
!   eselector         - user function to select open edges (OPTIONAL)
!
! OUTPUT
!   flow                 - the maximum flow from source to sink
!   position_mincutlabel - (OPTIONAL) partition the graph's vertices into two
!                          disjoint subsets that minimizes the total capacity
!                          of edges connecting the two subsets.
!                          "vertices/ipar" array item is labeled as
!                            0 - no flow through vertex (closed vertex)
!                            1 - source connected subset
!                            2 - sink connected subset
!                            3 - disconnected (open, but not accessible from
!                                the source vertex
!   position_flow        - (OPTIONAL)"edges/rpar" array item to save flow along
!                          the edge
!
! Output array items are updated for all graph objects.
!
      integer, parameter :: &
        CLOSED=0, SOURCE_REACHABLE=1, SINK_REACHABLE=2, DISCONNECTED=3
      integer, parameter :: NOT_DISCONNECTED=-1
      real(dp), allocatable :: forward_capacity(:), backward_capacity(:)
      integer, allocatable :: prev_edge(:), pair_edge(:)
      integer :: source_id, sink_id
      real(dp) :: additional_flow
      type(stack_t) :: added_edges
      logical, allocatable :: vmask0(:), emask0(:)

      ! Set up working arrays
      allocate(forward_capacity(this%nedges), backward_capacity(this%nedges))
        ! Remaining capacity for forward and backward flow, backward_capacity
        ! is used for undirected graphs only. For directed graphs reverse
        ! edges are temporarily added to the graph.
      allocate(prev_edge(this%nvertices))
        ! Keep track to the incoming edge id.

      ! Select open edges and vertices
      call graph_build_selection_masks(this, vmask0, emask0, &
        vselector=vselector, eselector=eselector, vmask_provided=vmask, &
        emask_provided=emask)

      block
        ! Capacity of all open edges (and edges connecting open vertices) is set
        ! to their capacity given in "edges/rpar" array.
        ! Capacity of closed edges set to zero (to disabling them)
        integer :: i, ia, ib
        do i=1, this%nedges
          forward_capacity(i) = 0.0_dp
          backward_capacity(i) = 0.0_dp
          if (emask0(i)) then
            ! these checks are being done in "build_selection_masks" and can be
            ! removed after testing TODO
            ia = get_index_from_handle(this, this%edges(i)%src_handle)
            ib = get_index_from_handle(this, this%edges(i)%dst_handle)
            ! orphaned edge
            if (ia==MAP_NULL .or. ib==MAP_NULL) then
              error stop 'graph_maxflow - selected edge has missing end-point'
              cycle
            end if
            ! edge with closed end-points
            if (.not. (vmask0(ia) .and. vmask0(ib))) then
              error stop 'graph_maxflow - selected edge has closed end-points'
              cycle
            end if

            ! edge is open and both end-points are also open
            forward_capacity(i) = this%edges(i)%rpar(position_capacity)
            backward_capacity(i) = forward_capacity(i)
          end if
        end do

        ! Verify sink and source vertices exist and are open
        source_id = get_index_from_handle(this, source)
        sink_id = get_index_from_handle(this, sink)
        if (source_id==MAP_NULL .or. sink_id==MAP_NULL) then
          error stop 'graph_max_flow - source/sink not found in graph'
        else if (source%handle_type/=VERTEX_HANDLE_TYPE .or. sink%handle_type/=VERTEX_HANDLE_TYPE) then
          error stop 'graph_max_flow - source/sink handles of unexpected type'
        else if (.not. vmask0(source_id)) then
          error stop 'graph_max_flow - source is not open'
        else if (.not. vmask0(sink_id)) then
          error stop 'graph_max_flow - sink is not open'
        else
          ! all assertions are ok
          continue
        end if
      end block

      ! Stack to store handles to temporarily added reverse edges.
      ! Used for directed graphs only.
      block
        type(handle_t) :: edge
        call added_edges%initialize(chunksize=size(transfer(edge,INTEGER_MOLD)))
      end block

      ! For directed graphs, reverse edges are added to the graph. This means
      ! that "forward_capacity" will be reallocated, "backward_capacity" will no
      ! longer be needed. Reference "pair_edge" will be used instead.
      if (this%is_directed_graph) then
        block
          type(handle_t) :: edge
          integer :: nreverse_edges, i, ireverse
          real(dp), allocatable :: tmp_forward_capacity(:)

          ! Count edges with non-zero capacity
          nreverse_edges = count(forward_capacity > 0.0_dp)

          ! For each non-zero capacity edge, a reverse edge is added
          allocate(tmp_forward_capacity(this%nedges+nreverse_edges), source=0.0_dp)
          allocate(pair_edge(this%nedges+nreverse_edges))

          ! Initialise pair_edge with self-pairs. Edges without an explicitly
          ! added residual reverse edge keep this mapping.
          do i=1,this%nedges
            pair_edge(i) = i
          end do

          ! Add temporary reverse edges for edges with non-zero capacity.
          do i=1,this%nedges
            if (.not. (forward_capacity(i)>0.0_dp)) cycle
            associate(e=>this%edges(i))
              edge = this%add_edge(e%dst_handle, e%src_handle, e%ipar, e%rpar)
            end associate
            ireverse = get_index_from_handle(this, edge)
            call added_edges%push(transfer(edge,INTEGER_MOLD))
            pair_edge(i) = ireverse
            pair_edge(ireverse) = i
            ! The capacity of reverse edges is initially set to zero as
            ! required by the algorithm.
            tmp_forward_capacity(i) = forward_capacity(i)
            tmp_forward_capacity(ireverse) = 0.0_dp
          end do
          call move_alloc(tmp_forward_capacity, forward_capacity)
          deallocate(backward_capacity)
          allocate(backward_capacity(0)) ! assert not be used later accidentaly
        end block
      else
        allocate(pair_edge(0)) ! array not needed for undirected graphs
      end if

      ! Make a complete BFS traversal to identify disocnnected vertices
      if (present(position_mincutlabel)) then
        block
          integer :: i
          call bfs_residual_search(this, forward_capacity, backward_capacity, &
              source_id, 0, prev_edge)
          do i=1,this%nvertices
            associate(label=>this%vertices(i)%ipar(position_mincutlabel))
              ! the labels are just temporary, will be relabeled later
              if (prev_edge(i)==MAP_NULL .and. i/=source_id) then
                ! this vertex could not be reached from source
                label = DISCONNECTED
              else
                label = NOT_DISCONNECTED
              end if
            end associate
          end do
        end block
      end if

      ! Initialize flow along edges (if required by user)
      if (present(position_flow)) then
        this%edges(1:this%nedges)%rpar(position_flow) = 0.0_dp
      end if

      ! The main loop of Edmonds-Karp
      ! Augment flow as long as path with non-zero capacity exists
      flow = 0.0_dp
      do
        ! Find shortest path using edges with non-zero remaining capacity
        call bfs_residual_search(this, forward_capacity, backward_capacity, &
            source_id, sink_id, prev_edge)
        if (prev_edge(sink_id)==MAP_NULL) exit
        ! The flow can be augmented. How much flow can we send?
        call process_path(this, forward_capacity, backward_capacity, &
            source_id, sink_id, prev_edge, pair_edge, additional_flow, .false.)
print *, 'Current flow is ', flow,'. Augmenting by ',additional_flow,'.'
        ! Update capacity of the network
        call process_path(this, forward_capacity, backward_capacity, &
            source_id, sink_id, prev_edge, pair_edge, additional_flow, .true., &
            position_flow)
        flow = flow + additional_flow
      end do

      ! Make minimum cut partition (if required by user)
      if (present(position_mincutlabel)) then
        block
          integer :: i

          ! unlimited traversal from the source
          call bfs_residual_search(this, forward_capacity, backward_capacity, &
              source_id, 0, prev_edge)

          do i=1,this%nvertices
            associate(label=>this%vertices(i)%ipar(position_mincutlabel))
              if (.not. vmask0(i)) then
                label = CLOSED
              else if (i==source_id) then
                label = SOURCE_REACHABLE
              else if (prev_edge(i)/=MAP_NULL) then
                label = SOURCE_REACHABLE
              else if (LABEL==DISCONNECTED) then
                ! open, could not be reached from source initially
                ! keep this label
                continue
              else
                ! open, reachable from the source initially, but unreachable
                ! in the residual network
                label = SINK_REACHABLE
              end if
            end associate
          end do
        end block
      end if

      ! Remove reverse edges added for directed graph
      block
        type(handle_t) :: edge
        do while(.not. added_edges%empty())
          edge = transfer(added_edges%pop(), edge)
          call this%remove_edge(edge)
        end do
      end block

    end subroutine graph_maxflow


    subroutine bfs_residual_search(this, forward_capacity, backward_capacity, &
        source_id, target_id, prev_edge)
      class(graph_t), intent(in) :: this
      real(dp), intent(in) :: forward_capacity(:), backward_capacity(:)
      integer, intent(in) :: source_id, target_id
      integer, intent(out) :: prev_edge(:)
!
! Breadth-first search of a residual network.
!
! Starting from "source_id", traverse edges with positive residual capacity
! and store the predecessor edge of each visited vertex in "prev_edge".
! If "target_id" is a valid vertex index, traversal stops after the target is
! reached. If "target_id" is non-positive, the complete reachable component
! is explored.
!
! The routine is used both for Edmonds-Karp augmenting path search and
! residual graph reachability analysis.
!
! Residual edge traversal:
!
!   Undirected graphs:
!     Traversing SRC -> DST uses forward_capacity.
!     Traversing DST -> SRC uses backward_capacity.
!
!   Directed graphs:
!     backward_capacity is not used.
!
! INPUT
!   this             - graph structure
!   forward_capacity - residual capacity in the forward direction
!   backward_capacity- residual capacity in the backward direction
!   source_id        - index of the starting vertex
!   target_id        - optional stopping vertex; non-positive value means
!                      unrestricted traversal
! OUTPUT
!   prev_edge        - predecessor edge used to reach each visited vertex;
!                      MAP_NULL for unvisited vertices and the source vertex
!
      type(queue_t) :: q
      integer :: current_id, iedge, ngb_id
      type(iterator_t) :: iterator

      call q%initialize(chunksize=size(transfer(current_id,INTEGER_MOLD)))
      call q%enqueue(transfer(source_id,INTEGER_MOLD))
      prev_edge = MAP_NULL

      do while(.not. q%empty())
        if (target_id > 0) then
          if (prev_edge(target_id)/=MAP_NULL) exit
        end if
        current_id = transfer(q%dequeue(), current_id)
        iterator = iterator_t()
        NGBS_LOOP: do while (this%vertices(current_id)%ngbs%has_next(iterator))
          call this%vertices(current_id)%ngbs%next(iterator, iedge)
          ngb_id = other_vertex_id(this, iedge, current_id)

          ! Skip edges with zero capacity
          if (ngb_id == get_index_from_handle(this,this%edges(iedge)%dst_handle)) then
            ! forward edge
            if (forward_capacity(iedge)<=0.0_dp) cycle
          else if (ngb_id == get_index_from_handle(this,this%edges(iedge)%src_handle)) then
            ! backward edge
            ! no backward edge can appear in directed graph
            if (this%is_directed_graph) error stop &
                'bfs_shortest_path - assertion for directed graph fails'
            if (backward_capacity(iedge)<=0.0_dp) cycle
          else
            error stop 'bfs_shortest_path - should not reach this branch'
          end if

          ! Skip edges going back to already traversed vertices
          if (ngb_id==source_id .or. prev_edge(ngb_id)/=MAP_NULL) cycle

          ! Add next node to the queue, mark which edge was used to come-in
          prev_edge(ngb_id) = iedge
          call q%enqueue(transfer(ngb_id,INTEGER_MOLD))
        end do NGBS_LOOP
      end do
      ! Now it is possible use "prev_edge(target_id)" to see if path from
      ! source to target exists and back-track the path back to source.
      ! The source vertex remains MAP_NULL as it has no predecessor.
    end subroutine bfs_residual_search


    subroutine process_path(this, forward_capacity, backward_capacity, &
        source_id, sink_id, prev_edge, pair_edge, additional_flow, &
        updating_flow, position_flow)
      class(graph_t), intent(inout) :: this
      real(dp), intent(inout) :: forward_capacity(:), backward_capacity(:)
      integer, intent(in) :: source_id, sink_id, prev_edge(:), pair_edge(:)
      real(dp), intent(inout) :: additional_flow
      logical, intent(in) :: updating_flow
      integer, intent(in), optional :: position_flow
!
! Back-track the path from sink to source and:
!  - find the bottleneck remaining capacity if "updating_flow==.false.",
!    or
!  - update remaining capacity along the path if "updating_flow==.true.".
!
      real(dp) :: capacity
      integer :: current_id, next_id
      logical :: is_forward_edge

      if (present(position_flow) .and. .not. updating_flow) &
          error stop 'position_flow argument can be given in update mode only'

      if (.not. updating_flow) additional_flow = huge(additional_flow)

      current_id = sink_id
      do while (prev_edge(current_id) /= MAP_NULL)
        next_id = other_vertex_id(this, prev_edge(current_id), current_id)

        if (next_id == get_index_from_handle(this,this%edges(prev_edge(current_id))%src_handle)) then
          is_forward_edge = .true.
          capacity = forward_capacity(prev_edge(current_id))
        else if (next_id == get_index_from_handle(this,this%edges(prev_edge(current_id))%dst_handle)) then
          ! backward edge
          ! no backward edge can appear in directed graph
          if (this%is_directed_graph) error stop &
              'process_path - assertion for directed graph fails'
          is_forward_edge = .false.
          capacity = backward_capacity(prev_edge(current_id))
        else
          error stop 'process_path - should not reach this branch'
        end if

        if (updating_flow) then
          ! Update capaciry mode
          if (this%is_directed_graph) then
            associate (fcap=>forward_capacity(prev_edge(current_id)), &
                bcap=>forward_capacity( pair_edge(prev_edge(current_id)) ) )
              fcap = fcap - additional_flow
              bcap = bcap + additional_flow
            end associate
            if (present(position_flow)) then
              associate (f=>this%edges(prev_edge(current_id))%rpar(position_flow), &
                  b=>this%edges(pair_edge(prev_edge(current_id)))%rpar(position_flow))
                f = f + additional_flow
                b = b - additional_flow
              end associate
            end if
          else
            ! undirected graph
            associate (fcap=>forward_capacity(prev_edge(current_id)), &
                bcap=>backward_capacity(prev_edge(current_id)))
              if (is_forward_edge) then
                fcap = fcap - additional_flow
                bcap = bcap + additional_flow
              else
                fcap = fcap + additional_flow
                bcap = bcap - additional_flow
              end if
            end associate
            if (present(position_flow)) then
              associate (f=>this%edges(prev_edge(current_id))%rpar(position_flow))
                if (is_forward_edge) then
                  f = f + additional_flow
                else
                  f = f - additional_flow
                end if
              end associate
            end if
          end if
        else
          ! Looking for the bottleneck mode
          if (capacity < additional_flow) additional_flow = capacity
        end if

        current_id = next_id
      end do

      ! verify source reached
      if (current_id /= source_id) error stop &
          'process_path - could not reach source'
    end subroutine process_path


    subroutine graph_maxflow_multiple(this, sources, sinks, &
        position_capacity, flow, position_mincutlabel, position_flow, &
        vmask, emask, vselector, eselector)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: sources(:), sinks(:)
      integer, intent(in) :: position_capacity
      real(dp), intent(out) :: flow
      integer, intent(in), optional :: position_mincutlabel
      integer, intent(in), optional :: position_flow
      logical, intent(in), optional :: vmask(:), emask(:)
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
!
! Maximum flow using multiple sources and sinks
!
      type(handle_t) :: super_source, super_sink, edge
      type(stack_t) :: added_edges
      integer :: nvertices0, nedges0
      real(dp) :: total_capacity
      logical, allocatable :: vmask0(:), emask0(:)

      ! Verify the source and sink lists:
      ! - at least one source and one sink are present
      ! - all handles are of VERTEX_HANDLE_TYPE
      ! - all handles are valid and unique
      block
        integer, allocatable :: listed_count(:)
        integer :: i, iv

        if (size(sources)<1 .or. size(sinks)<1) error stop &
            'graph_maxflow_multiple - zero source/sink vertices'

        if (any(sources%handle_type /= VERTEX_HANDLE_TYPE) .or. &
            any(sinks%handle_type /= VERTEX_HANDLE_TYPE)) error stop &
            'graph_maxflow_multiple - all source and sink handles must be vertices'

        allocate(listed_count(this%nvertices), source=0)
        do i=1, size(sources)
          iv = get_index_from_handle(this, sources(i))
          if (iv==MAP_NULL) error stop &
              'graph_maxflow_multiple - a source handle not found in graph'
          listed_count(iv) = listed_count(iv)+1
        end do
        do i=1, size(sinks)
          iv = get_index_from_handle(this, sinks(i))
          if (iv==MAP_NULL) error stop &
              'graph_maxflow_multiple - a sink handle not found in graph'
          listed_count(iv) = listed_count(iv)+1
        end do
        if (any(listed_count>1)) error stop &
          'graph_maxflow_multiple - source/sink verticies must be unique'
      end block

      ! Save number of objects for assertion at the end
      nvertices0 = this%nvertices
      nedges0 = this%nedges

      ! Select open edges and vertices
      call graph_build_selection_masks(this, vmask0, emask0, &
          vselector=vselector, eselector=eselector, &
          vmask_provided=vmask, emask_provided=emask)

      ! Sum the capacity over all open edges to be used as the
      ! capacity of added edges connecting super nodes.
      total_capacity = &
          sum(this%edges(1:this%nedges)%rpar(position_capacity),mask=emask0)

      ! Add super-source and super-sink and connect them to sources and sinks.
      block
        integer :: v_ipar(VSIZE_IPAR), e_ipar(ESIZE_IPAR), i
        integer :: nopen_sources, nopen_sinks
        real(dp) :: v_rpar(VSIZE_RPAR), e_rpar(ESIZE_RPAR)

        v_ipar = 0
        v_rpar = 0.0_dp
        e_ipar = 0
        e_rpar = 0.0_dp
        e_rpar(position_capacity) = total_capacity
        super_source = this%add_vertex(v_ipar, v_rpar)
        super_sink = this%add_vertex(v_ipar, v_rpar)
        call added_edges%initialize(chunksize=size(transfer(edge,INTEGER_MOLD)))
        nopen_sources = 0
        do i=1, size(sources)
          ! if source is closed, do not add the connection
          if (.not. vmask0(get_index_from_handle(this, sources(i)))) cycle
          edge = this%add_edge(super_source, sources(i), e_ipar, e_rpar)
          call added_edges%push(transfer(edge,INTEGER_MOLD))
          nopen_sources = nopen_sources+1
        end do
        nopen_sinks = 0
        do i=1, size(sinks)
          ! if sink is closed, do not add the connection
          if (.not. vmask0(get_index_from_handle(this, sinks(i)))) cycle
          edge = this%add_edge(sinks(i), super_sink, e_ipar, e_rpar)
          call added_edges%push(transfer(edge,INTEGER_MOLD))
          nopen_sinks = nopen_sinks+1
        end do

        if (nopen_sinks==0 .or. nopen_sources==0) &
          print '("maxflow_multiple WARNING - zero flow as all sources or sinks closed")'
      end block

      ! Extend masks to include super-source, super-sink and
      ! their connecting edges.
      block
        logical, allocatable :: vmask_tmp(:), emask_tmp(:)
        allocate(vmask_tmp(size(vmask0)+2), source=.true.)
        allocate(emask_tmp(size(emask0)+added_edges%size()), source=.true.)
        vmask_tmp(1:size(vmask0)) = vmask0
        emask_tmp(1:size(emask0)) = emask0
        call move_alloc(vmask_tmp, vmask0)
        call move_alloc(emask_tmp, emask0)
      end block

      ! Max-flow
      call graph_maxflow( &
          this, super_source, super_sink, position_capacity, flow, &
          position_mincutlabel=position_mincutlabel, &
          position_flow=position_flow, &
          vmask=vmask0, emask=emask0)

      ! Remove added edges/vertices and assert number of objects did not change
      do while (.not. added_edges%empty())
        call this%remove_edge(transfer(added_edges%pop(),edge))
      end do
      call this%remove_vertex(super_sink)
      call this%remove_vertex(super_source)
      if (this%nvertices/=nvertices0) error stop &
          'graph_maxflow_multiple - number of vertices changed (internal error)'
      if (this%nedges/=nedges0) error stop &
          'graph_maxflow_multiple - number of edges changed (internal error)'

    end subroutine graph_maxflow_multiple


    ! -------------------------------
    ! Betweenness (Brandes algorithm)
    ! -------------------------------
    subroutine graph_betweenness(this, position_cost, position_eb, position_vb, &
        is_normalized, vselector, eselector, vmask, emask)
      class(graph_t), intent(inout) :: this
      integer, intent(in), optional :: position_cost
      integer, intent(in), optional :: position_eb, position_vb
      logical, intent(in), optional :: is_normalized
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      logical, intent(in), optional :: vmask(:), emask(:)
!
!  Compute edge and vertex betweenness centrality using Brandes' algorithm.
!  Shortest paths are computed using Dijkstra's algorithm and edge costs stored
!  in "edges(:)%rpar(position_cost)".
!  Optional vertex and edge selectors restrict the computation to a selected
!  subgraph. Betweenness values are accumulated only for selected vertices and
!  edges. Vertices and edges excluded from the selection are left unchanged.
!
!  For undirected graphs, the accumulated scores are divided by two because every
!  shortest path is encountered twice.
!
!  If requested, the scores are normalized according to the number of selected
!  vertices.
!
      logical, allocatable :: vmask0(:), emask0(:)
      integer :: i, id_s
      real(dp), allocatable :: delta(:)
        ! dependency of the source on vertex "v"
      type(stack_t), allocatable :: prev(:)
        ! list of immediate predecessors of "v" on shortest paths
      integer(I8B), allocatable :: sigma(:)
        ! the number of unique shortest paths from "s" to "v"
      type(stack_t) :: stack

      ! Dijkstra search working storage
      type(pqueue_t) :: pqueue_dijkstra
      real(dp), allocatable :: dist_dijkstra(:)
        ! shortest distance from "s" to "v" in weighted graphs
      type(pqueue_handle_t), allocatable :: phas(:)
        ! handles to vertices added to Dijksta's priority queue
      logical, allocatable :: visited(:)
        ! denote visited nodes during Dijksta search

      ! BFS working storage
      integer, allocatable :: dist_bfs(:)
        ! shortest distance from "s" to "v" in unweighted graphs
      type(queue_t) :: queue_bfs

      ! Nothing to do if both "position_eb" and "position_vb" are ommitted
      if (.not. (present(position_eb) .or. present(position_vb))) return

      ! Mark selected edges and count selected_vertices
      call graph_build_selection_masks(this, vmask0, emask0, &
          vselector=vselector, eselector=eselector, &
          vmask_provided=vmask, emask_provided=emask)

      ! Initialize storage for Dijkstra/BFS
      allocate(delta(this%nvertices))
      allocate(sigma(this%nvertices))
      allocate(prev(this%nvertices))
      do i=1, this%nvertices
        call prev(i)%initialize(size(transfer(id_s,INTEGER_MOLD)))
      end do
      call stack%initialize(chunksize=size(transfer(id_s,INTEGER_MOLD)))

      if (present(position_cost)) then
        ! Weighted graph: initialize storage for Dijkstra search
        allocate(dist_dijkstra(this%nvertices))
        allocate(phas(this%nvertices), visited(this%nvertices))
        call pqueue_dijkstra%initialize(chunksize=size(transfer(id_s,INTEGER_MOLD)))
      else
        ! Unweighted graph: initialize storage for BFS
        allocate(dist_bfs(this%nvertices))
        call queue_bfs%initialize(chunksize=size(transfer(id_s,INTEGER_MOLD)))
      end if

      ! Initialize edge and vertex betweeness
      if (present(position_eb)) then
        where (emask0) &
            this%edges(1:this%nedges)%rpar(position_eb) = 0.0_dp
      end if
      if (present(position_vb)) then
        where (vmask0) &
            this%vertices(1:this%nvertices)%rpar(position_vb) = 0.0_dp
      end if

      ! Main loop over all source vectors
      SRC_LOOP: do id_s=1, this%nvertices
if (mod(id_s,5000)==0) print '("Source is ",i0," out of ",i0)', id_s, this%nvertices
        if (.not. vmask0(id_s)) cycle

        ! STEP 1 - Shortest paths search
        if (present(position_cost)) then
          call betweenness_dijkstra(this, id_s, emask0, position_cost, sigma, &
              stack, prev, pqueue_dijkstra, dist_dijkstra, visited, phas)
        else
          call betweenness_bfs(this, id_s, emask0, sigma, stack, prev, queue_bfs, dist_bfs)
        end if

        ! STEP 2 - Backward pass
        delta = 0.0_dp
        BACKPASS_LOOP: do while(.not. stack%empty())
          block
            integer :: id_v, id_u, iedge
            real(dp) :: delta_edge
            id_v = transfer(stack%pop(), id_v)
            do while (.not. prev(id_v)%empty())
              iedge = transfer(prev(id_v)%pop(), iedge)
              id_u = other_vertex_id(this, iedge, id_v)
              delta_edge = real(sigma(id_u),dp)/real(sigma(id_v),dp)*(1.0_dp+delta(id_v))
              ! accumulate global edge beteenness
              if (present(position_eb)) then
                associate(c=>this%edges(iedge)%rpar(position_eb))
                  c = c + delta_edge
                end associate
              end if
              ! pass it to a predecessor node
              delta(id_u) = delta(id_u) + delta_edge
            end do
            ! accumulate vertex betweenness
            if (id_v /= id_s .and. present(position_vb)) then
              associate(c=>this%vertices(id_v)%rpar(position_vb))
                c = c + delta(id_v)
              end associate
            end if
          end block
        end do BACKPASS_LOOP

      end do SRC_LOOP

      ! Divide by two for undirected graphs (all paths were counted twice)
      if (.not. this%is_directed_graph) then
        if (present(position_eb)) then
          associate(eb=>this%edges(1:this%nedges)%rpar(position_eb))
            where (emask0) eb = 0.5_dp * eb
          end associate
        end if
        if (present(position_vb)) then
          associate(vb=>this%vertices(1:this%nvertices)%rpar(position_vb))
            where (vmask0) vb = 0.5_dp * vb
          end associate
        end if
      end if

      ! Normalize the scores (if asked for by an user)
      block
        logical :: is_normalized0
        real(dp) :: v_denominator, e_denominator
        integer :: n
        is_normalized0 = .false. ! default behaviour
        if (present(is_normalized)) is_normalized0 = is_normalized
        if (is_normalized0) then
          n = count(vmask0)
          if (n>=3 .and. present(position_vb)) then
            ! Vertices
            ! directed graphs: total ordered pairs of nodes (paths s-->t  and
            !                  t-->s are counted as separate pairs).
            ! undirected graphs: total unordered pairs of nodes in the graph,
            !                    excluding the target node itself.
            v_denominator = real((n-1)*(n-2),dp)
            if (.not. this%is_directed_graph) v_denominator = v_denominator * 0.5_dp
            associate(vb=>this%vertices(1:this%nvertices)%rpar(position_vb))
              where (vmask0) vb = vb / v_denominator
            end associate
          end if
          if (n>=2 .and. present(position_eb)) then
            ! Edges
            ! directed graphs: total possible unique ordered pairs of nodes in
            !                  the entire graph.
            ! undirected graphs: total possible unique unordered pairs of nodes
            !                    in the entire graph.
            e_denominator = real(n*(n-1),dp)
            if (.not. this%is_directed_graph) e_denominator = e_denominator * 0.5_dp
            associate(eb=>this%edges(1:this%nedges)%rpar(position_eb))
              where(emask0) eb = eb / e_denominator
            end associate
          end if
        end if
      end block

      ! Clean-up
      ! got run time error (compiler bug?), explicit deallocation solved this
      if (present(position_cost)) deallocate(phas)

    end subroutine graph_betweenness


    subroutine betweenness_dijkstra(this, id_s, emask0, position_cost, sigma, stack, prev, pqueue, dist, visited, phas)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: id_s, position_cost
      logical, intent(in) :: emask0(:)
      type(pqueue_t), intent(inout) :: pqueue
      type(stack_t), intent(inout) :: stack
      real(dp), intent(out) :: dist(:)
      type(stack_t), intent(inout) :: prev(:)
      integer(I8B), intent(out) :: sigma(:)
      logical, intent(out) :: visited(:)
      type(pqueue_handle_t), intent(inout) :: phas(:)
!
! Shortest paths search in weighted graph
!
      integer :: id_v, id_u, edge_uv
      real(dp) :: dist_to_v
      type(iterator_t) :: iterator
      integer, parameter :: DIST_SHORTER=1, DIST_SAME=0, DIST_LONGER=-1

      ! Initialize Dijkstra's search from source
      dist = huge(dist)
      do id_u=1, this%nvertices
        call prev(id_u)%clear()
        phas(id_u) = pqueue_handle_t()
      end do
      sigma = 0_I8B
      visited = .false.

      dist(id_s) = 0.0_dp
      sigma(id_s) = 1_I8B
      phas(id_s) = pqueue%insert(transfer(id_s,INTEGER_MOLD), dist(id_s))

      DJIKSTRA_LOOP: do while (.not. pqueue%empty())
        ! dequeue vertex and push it to the stack for later use
        id_u = transfer(pqueue%pop(), id_u)
        call stack%push(transfer(id_u,INTEGER_MOLD))
        visited(id_u) = .true.

        iterator = iterator_t()
        NGB_LOOP: do while(this%vertices(id_u)%ngbs%has_next(iterator))
          call this%vertices(id_u)%ngbs%next(iterator, edge_uv)
          if (.not. emask0(edge_uv)) cycle
          id_v = other_vertex_id(this, edge_uv, id_u)
          if (visited(id_v)) cycle

          dist_to_v = dist(id_u) + this%edges(edge_uv)%rpar(position_cost)
          select case(compare_dist(dist_to_v, dist(id_v)))
          case(DIST_SHORTER)
            dist(id_v) = dist_to_v
            call prev(id_v)%clear()
            call prev(id_v)%push(transfer(edge_uv,INTEGER_MOLD))
            sigma(id_v) = sigma(id_u)
            if (pqueue%contains(phas(id_v))) then
              call pqueue%update_priority(phas(id_v), dist_to_v)
            else
              phas(id_v) = pqueue%insert(transfer(id_v,INTEGER_MOLD), dist_to_v)
            end if
          case(DIST_SAME)
            call prev(id_v)%push(transfer(edge_uv,INTEGER_MOLD))
            sigma(id_v) = sigma(id_v) + sigma(id_u)
          end select
        end do NGB_LOOP
      end do DJIKSTRA_LOOP

    contains
      integer function compare_dist(new, old)
        real(dp), intent(in) :: old, new
        real(dp), parameter :: REL_TOL = 1.0e5_dp * epsilon(1.0_dp)
        real(dp) :: tol
        tol = REL_TOL * max(1.0_dp, abs(new), abs(old))
        if (abs(old-new) < tol) then
          compare_dist = DIST_SAME
        else if (new < old) then
          compare_dist = DIST_SHORTER
        else
          compare_dist = DIST_LONGER
        end if
      end function
    end subroutine betweenness_dijkstra


    subroutine betweenness_bfs(this, id_s, emask0, sigma, stack, prev, queue, dist)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: id_s
      logical, intent(in) :: emask0(:)
      integer(I8B), intent(out) :: sigma(:)
      type(stack_t), intent(inout) :: stack
      type(stack_t), intent(inout) :: prev(:)
      type(queue_t), intent(inout) :: queue
      integer, intent(out) :: dist(:)
!
! Shortest paths search in unweighted graphs.
!
      integer :: id_u, id_v, edge_uv
      type(iterator_t) :: iterator

      ! initialize
      dist = -1
      sigma = 0_I8B
      do id_u = 1, this%nvertices
        call prev(id_u)%clear()
      end do

      dist(id_s) = 0
      sigma(id_s) = 1_I8B
      call queue%enqueue(transfer(id_s,INTEGER_MOLD))

      BFS_LOOP: do while(.not. queue%empty())
        ! dequeue vertex and push it to the stack for later use
        id_u = transfer(queue%dequeue(), id_u)
        call stack%push(transfer(id_u,INTEGER_MOLD))

        iterator = iterator_t()
        NGB_LOOP: do while(this%vertices(id_u)%ngbs%has_next(iterator))
          call this%vertices(id_u)%ngbs%next(iterator, edge_uv)
          if (.not. emask0(edge_uv)) cycle
          id_v = other_vertex_id(this, edge_uv, id_u)
          if (dist(id_v)<0) then ! found for the first time
            call queue%enqueue(transfer(id_v,INTEGER_MOLD))
            dist(id_v) = dist(id_u)+1
          end if
          if (dist(id_v) == dist(id_u)+1) then ! shortest path to "v" via "u"
            sigma(id_v) = sigma(id_v)+sigma(id_u)
            call prev(id_v)%push(transfer(edge_uv,INTEGER_MOLD))
          end if
        end do NGB_LOOP
      end do BFS_LOOP

    end subroutine betweenness_bfs


    ! --------------------------------------------
    ! Stoer-Wagner algorithm for min-cut partition
    ! --------------------------------------------
    subroutine graph_mincut(this, position_weight, mincut, s_list, t_list, &
        vmask, emask, vselector, eselector)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: position_weight
      real(dp), intent(out) :: mincut
      type(handle_t), intent(out), allocatable :: s_list(:), t_list(:)
      logical, intent(in), optional :: vmask(:), emask(:)
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
!
! Find min-cut
!
      type(graph_t) :: g
      type(handle_t) :: handle
      type(stack_t), allocatable :: original_vertices(:)
      type(stack_t) :: mincut_vertices
      integer, allocatable :: original_labels(:)
      type(pqueue_handle_t), allocatable :: pqueue_handles(:)
      type(pqueue_t) :: scoreboard
      integer, parameter :: LAB_NOT_SELECTED=0, LAB_SET_S=1, LAB_SET_T=2

      ! Stoer-Wagner algorithm works for undirected graphs only
      if (this%is_directed_graph) error stop &
          'graph_mincut - algorithm requires an undirected graph only'

      ! Prepare the working graph and the array of stacks with handles to the
      ! vertices in the original graph.
      block
        integer :: i, k
        type(handle_t), allocatable :: vertices0(:)

        allocate(original_labels(this%nvertices), source=LAB_NOT_SELECTED)
        call this%copy(g, vselector=vselector, eselector=eselector, &
            vmask=vmask, emask=emask, new_vertices=vertices0)
        if (g%nvertices < 2) error stop &
            'graph_mincut - at least two vertices are required'
        allocate(original_vertices(g%nvertices))
        do i=1, g%nvertices
          if (g%vertices(i)%handle%index_to_map /= i) error stop &
            'graph_mincut - internal assertion fails (1)'
          call original_vertices(i)%initialize( &
              chunksize=size(transfer(handle,INTEGER_MOLD)))
        end do
        i = 1
        do k=1, size(vertices0)
          if (get_index_from_handle(this, vertices0(k))==MAP_NULL) cycle
          ! The k-th vertex in the original graph is now i-th vertex in the
          ! working graph. Push the handle the vertex is known by in the
          ! original graph to the stack maped with the working graph.
          call original_vertices(i)%push(transfer(vertices0(k),INTEGER_MOLD))
          original_labels(k) = LAB_SET_S
          i = i+1
        end do
        if (i-1/=g%nvertices) error stop &
            'graph_mincut - internal assertion fails (2)'
      end block

      ! The priority queue stores the handles of vertices not yet moved
      ! to the set A, the vertex connectivity with vertices already present
      ! in the set A is stored as the priority in the queue.
      call scoreboard%initialize( &
          chunksize=size(transfer(handle,INTEGER_MOLD)), ordering=PQUEUE_MAX)
      allocate(pqueue_handles(g%nvertices))

      ! Initialize global mincut value and corresponding vertices
      call mincut_vertices%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
      mincut = huge(mincut)

      block
        type(handle_t) :: s_handle, t_handle
        real(dp) :: mincut_now
        MAIN_LOOP: do
          ! phase 1
          call find_st(g, scoreboard, position_weight, s_handle, pqueue_handles)
          t_handle = transfer(scoreboard%pop(top_priority=mincut_now), t_handle)
          if (.not. scoreboard%empty()) error stop &
              'graph_mincut - internal assertion fails (3)'

          ! update min-cut value and vertices list if a lower value found
          if (mincut_now < mincut) then
            mincut = mincut_now
            mincut_vertices = original_vertices(t_handle%index_to_map)
          end if

          ! phase 2 - contract vertices S and T
          if (g%nvertices <= 2) exit MAIN_LOOP
          call contract_st(g, position_weight, s_handle, t_handle, original_vertices)
        end do MAIN_LOOP
      end block

      ! Consume mincut_vertices to label vertices from the winning mincut
      do while (.not. mincut_vertices%empty())
        handle = transfer(mincut_vertices%pop(), handle)
        associate(lab=>original_labels(get_index_from_handle(this, handle)))
          if (lab /= LAB_SET_S) error stop &
              'graph_mincut - internal assertion fails (4)'
          lab = LAB_SET_T
        end associate
      end do

      ! Make s_list and t_list
      block
        integer :: is, it, k
        allocate(s_list(count(original_labels==LAB_SET_S)))
        allocate(t_list(count(original_labels==LAB_SET_T)))
        is = 1
        it = 1
        do k=1, size(original_labels)
          if (original_labels(k)==LAB_SET_S) then
            s_list(is) = this%vertices(k)%handle
            is = is+1
          else if (original_labels(k)==LAB_SET_T) then
            t_list(it) = this%vertices(k)%handle
            it = it+1
          end if
        end do
        if (is-1/=size(s_list) .or. it-1/=size(t_list)) error stop &
            'graph_mincut - internal assertion fails (5)'
      end block

      ! Clean-up (auto-deallocation did not work?)
      deallocate(pqueue_handles)

print '("graph_mincut: Min weight is ",g0,".",/,"&
    &Set S contains ",i0," nodes, set T contains ",i0," nodes and ",i0," nodes were not selected")', &
    mincut, size(t_list), size(s_list), count(original_labels==LAB_NOT_SELECTED) 

    end subroutine graph_mincut


    subroutine find_st(g, scoreboard, position_weight, s_handle, phandles)
      type(graph_t), intent(in) :: g
      type(pqueue_t), intent(inout) :: scoreboard
      integer, intent(in) :: position_weight
      type(handle_t), intent(out) :: s_handle
      type(pqueue_handle_t), intent(inout) :: phandles(:)
!
! Move vertices to set A as defined by Stoer-Wagner algorithm.
! The working priority queue contains vertices outside of set A.
!
      integer :: s_imap, s_id, w_imap, w_id, ie
      type(iterator_t) :: iterator

      ! Defensive checks
      if (.not. scoreboard%empty()) error stop &
          'find_st - the queue is not empty (internal error)'
      if (g%nvertices < 2) error stop &
          'find_st - less than two vertices in graph (internal error)'

      ! Store all vertices to the priority queue...
      do s_id=1, g%nvertices
        s_imap = g%vertices(s_id)%handle%index_to_map
        phandles(s_imap)=scoreboard%insert( &
            transfer(g%vertices(s_id)%handle,INTEGER_MOLD), priority=0.0_dp)
      end do

      ! ...and iteratively remove the vertex S with the highest connectivity
      do while(scoreboard%size() > 1)
        s_handle = transfer(scoreboard%pop(), s_handle)
        s_imap = s_handle%index_to_map
        s_id = get_index_from_handle(g, s_handle)

        ! as S now becomes part of set A, increase connectivity
        ! of all vertices W that are neighbors of S and are outside of set A 
        iterator = iterator_t()
        do while (g%vertices(s_id)%ngbs%has_next(iterator))
          call g%vertices(s_id)%ngbs%next(iterator, ie)
          w_id = other_vertex_id(g, ie, s_id)
          w_imap = g%vertices(w_id)%handle%index_to_map
          associate (handle=>phandles(w_imap))
            if (scoreboard%contains(handle)) &
                call scoreboard%update_priority(handle, &
                    new_priority = scoreboard%priority(handle) + &
                    g%edges(ie)%rpar(position_weight))
          end associate
        end do

        ! on return, the queue contains the last vertex T
        if (scoreboard%size()==1) exit
      end do

    end subroutine find_st


    subroutine contract_st(g, position_weight, s_handle, t_handle, voriginal)
      class(graph_t), intent(inout) :: g
      integer, intent(in) :: position_weight
      type(handle_t), intent(in) :: s_handle, t_handle
      type(stack_t), intent(inout) :: voriginal(:)
!
! Contract the graph by joining S and T vertices. From the data storage point
! of view, vertex S becomes the super-vertex ST, while vertex T is removed
! after its assets are transferred to vertex S.
!
      type(iterator_t) :: iterator
      type(handle_t) :: handle
      integer :: ie, isw, s_id, t_id, w_id, e_ipar(ESIZE_IPAR)
      real(dp) :: e_rpar(ESIZE_RPAR)

      ! The actual position of respective vertices in graph arrays
      t_id = get_index_from_handle(g, t_handle)
      s_id = get_index_from_handle(g, s_handle)

      ! Transfer all original vertices associated with vertex T to vertex S.
      associate (T=>voriginal(t_handle%index_to_map), &
                 S=>voriginal(s_handle%index_to_map) )
        do while (.not. T%empty())
          call S%push(T%pop())
        end do
      end associate

      ! Loop over all T's connections, remove them and transfer their weights
      ! to vertex S.
      iterator = iterator_t()
      do while (g%vertices(t_id)%ngbs%has_next(iterator))
        ! Obtain the index of the T--W edge. Because this edge will be removed,
        ! move then the iterator one item back.
        call g%vertices(t_id)%ngbs%next(iterator, ie)
        call g%vertices(t_id)%ngbs%back(iterator)
        ! Before edge is removed, get weight of removed edge and the index of
        ! neighbour W. The replacement edge, that may be added, will inherit the
        ! properties of the removed edge.
        w_id = other_vertex_id(g, ie, t_id)
        e_ipar = g%edges(ie)%ipar
        e_rpar = g%edges(ie)%rpar ! stores the removed edge's weight
        call g%remove_edge(g%edges(ie)%handle)

        ! If the edge is S--T, it is just removed. For other edges:
        ! - edge S--W with the weight of the removed edge T--W is added, or,
        ! - the weight of S--W edge is increased by the weight of the removed
        !   edge T--W if edge S--W already exists.
        if (w_id == s_id) then
          continue
        else
          isw = g%find_edge_id(s_id, w_id)
          if (isw == MAP_NULL) then
            ! add S--W edge as the replacement of the removed T--W edge
            handle = g%add_edge(s_handle, g%vertices(w_id)%handle, e_ipar, e_rpar)
          else
            ! increase capacity of the existing S--W edge
            associate(weight=>g%edges(isw)%rpar(position_weight))
              weight = weight + e_rpar(position_weight)
            end associate
          end if
        end if
      end do

      ! Vertex T should be isolated and can be removed
      if (g%vertices(t_id)%ngbs%size()/=0) error stop &
          'join_st - t still has connections (internal error)'
      call g%remove_vertex(t_handle)
    end subroutine contract_st


    ! ----------------------
    ! Edge / Vertex selector
    ! ----------------------
    subroutine graph_build_selection_masks(this, vmask, emask, vmask_provided, &
        emask_provided, vselector, eselector)
      class(graph_t), intent(in) :: this
      logical, intent(inout), allocatable :: vmask(:), emask(:)
      logical, intent(in), optional :: vmask_provided(:), emask_provided(:)
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
!
! Build vertex and edge selection masks from user-provided masks or selector
! functions.
!
! This routine creates logical masks identifying the vertices and edges to be
! considered by graph algorithms. Selection can be specified either by providing
! logical mask arrays or by supplying selector functions. If neither is provided,
! all vertices and edges are selected.
!
! A mask array and a selector function cannot be provided simultaneously for the
! same entity type. Edge selection is additionally restricted by vertex
! selection: an edge is selected only if both of its endpoint vertices are
! selected. (Every selected edge connects two selected endpoints)
!
! Existing allocated output masks are reused if they have the correct size;
! otherwise they are reallocated.
!
! Arguments:
!   this            - Graph object.
!   vmask, emask    - Output vertex and edge selection masks.
!   vmask_provided,
!   emask_provided  - Optional user-provided vertex and edge masks.
!   vselector,
!   eselector       - Optional vertex and edge selector functions.
!
      integer :: ia, ib, ie

      ! Only mask array or selector function can by provided, not both
      if (present(emask_provided) .and. present(eselector)) error stop &
          'build_selection_masks - both edge mask and selector function provided'
      if (present(vmask_provided) .and. present(vselector)) error stop &
          'build_selection_masks - both vertex mask and selector function provided'

      ! Assert provided masks have correct size
      if (present(vmask_provided)) then
        if (size(vmask_provided)/=this%nvertices) error stop &
            'build_selection_masks - vmask_provided has wrong size'
      end if
      if (present(emask_provided)) then
        if (size(emask_provided)/=this%nedges) error stop &
            'build_selection_masks - emask_provided has wrong size'
      end if

      ! Verify output mask arrays allocated to correct size.
      ! Allocate / reallocate if needed
      if (allocated(vmask)) then
        if (size(vmask) /= this%nvertices) deallocate(vmask)
      end if
      if (allocated(emask)) then
        if (size(emask) /= this%nedges) deallocate(emask)
      end if
      if (.not. allocated(vmask)) allocate(vmask(this%nvertices))
      if (.not. allocated(emask)) allocate(emask(this%nedges))

      ! Build selection array for verices first...
      if (present(vmask_provided)) then
        vmask = vmask_provided
      else if (present(vselector)) then
        do ia=1, this%nvertices
          vmask(ia) = vselector(this, this%vertices(ia))
        end do
      else
        vmask = .true.
      end if

      ! ...and build selection array for edges next
      if (present(emask_provided)) then
        emask = emask_provided
      else if (present(eselector)) then
        do ie=1, this%nedges
          emask(ie) = eselector(this, this%edges(ie))
        end do
      else
        emask = .true.
      end if

      do ie=1, this%nedges
        ! close edges with closed vertex as one of its ends
        if (.not. emask(ie)) cycle
        ia = get_index_from_handle(this, this%edges(ie)%src_handle)
        ib = get_index_from_handle(this, this%edges(ie)%dst_handle)
        if (ia==MAP_NULL .or. ib==MAP_NULL) then
          ! dangling edge not selected
          emask(ie) = .false.
        else if (.not. (vmask(ia) .and. vmask(ib))) then
          ! one of vertices is closed
          emask(ie) = .false.
        end if
      end do
    end subroutine graph_build_selection_masks


    function graph_select_vertices(this, vselector, vmask) result(handles)
      class(graph_t), intent(in) :: this
      procedure(is_vertex_selected), optional :: vselector
      logical, intent(in), optional :: vmask(:)
      type(handle_t), allocatable :: handles(:)
!
! Return handles to all vertices selected by vselector.
! If vselector is absent, handles to all vertices are returned.
! The order of handles follows the current internal vertex ordering.
!
      logical, allocatable :: vmask0(:)
      integer :: i, k

      if (present(vselector) .and. present(vmask)) error stop &
          'graph_select_vertices - both mask and selector function provided'

      allocate(vmask0(this%nvertices))
      if (present(vselector)) then
        do i=1, this%nvertices
          vmask0(i) = vselector(this, this%vertices(i))
        end do
      else if (present(vmask)) then
        if (size(vmask)/=this%nvertices) error stop &
            'graph_select_vertices - vmask has wrong size'
        vmask0 = vmask
      else
        vmask0 = .true.
      end if

      allocate(handles(count(vmask0)))
      k = 1
      do i=1, this%nvertices
        if (.not. vmask0(i)) cycle
        handles(k) = this%vertices(i)%handle
        k = k+1
      end do

      ! Defensive assertion
      if (k-1 /= size(handles)) error stop &
        'graph_select_vertices - internal assertion failed'
    end function graph_select_vertices


    function graph_select_edges(this, eselector, emask) result(handles)
      class(graph_t), intent(in) :: this
      procedure(is_edge_selected), optional :: eselector
      logical, intent(in), optional :: emask(:)
      type(handle_t), allocatable :: handles(:)
!
! Return handles to all edges selected by eselector or emask.
! If both eselector and emask are absent, handles to all edges are returned.
! The order of handles follows the current internal edge ordering.
!
      logical, allocatable :: emask0(:)
      integer :: i, k

      if (present(eselector) .and. present(emask)) error stop &
          'graph_select_edges - both mask and selector function provided'

      allocate(emask0(this%nedges))
      if (present(eselector)) then
        do i=1, this%nedges
          emask0(i) = eselector(this, this%edges(i))
        end do
      else if (present(emask)) then
        if (size(emask)/=this%nedges) error stop &
            'graph_select_edges - emask has wrong size'
        emask0 = emask
      else
        emask0 = .true.
      end if

      allocate(handles(count(emask0)))
      k = 1
      do i=1, this%nedges
        if (.not. emask0(i)) cycle
        handles(k) = this%edges(i)%handle
        k = k+1
      end do

      ! Defensive assertion
      if (k-1 /= size(handles)) error stop &
        'graph_select_edges - internal assertion failed'
    end function graph_select_edges

  end module graph_mod
