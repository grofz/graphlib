  module graph_mod
    use iso_fortran_env, only : DP => real64, I1B => int8, I8B => int64, output_unit
    use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR
    use conts_mod, only : queue_t, stack_t, pqueue_t, pqueue_handle_t=>handle_t, &
      PQUEUE_MIN, PQUEUE_MAX
    use graph_adjlist_mod, only : adjlist_t, iterator_t
    implicit none (type, external)
    private

    ! Workaround to gfortran bug (procedures used by submodules)
    public other_vertex_id, get_index_from_handle

    ! Parametrized derived type (PDT) not working reliably with compilers.
    ! To avoid PDT, array sizes required for the actual implementation
    ! are hardcoded in "graph_user.f90" and imported as ?SIZE_?PAR named
    ! constants
    !
    ! Alternativelly, array sizes can be hardcoded here
    ! integer, parameter :: VSIZE_IPAR=?, VSIZE_RPAR=?, ESIZE_IPAR=?, ESIZE_RPAR=?

    ! Exported constants used by module procedures (graph_mincut)
    integer, parameter, public :: &
        MINCUT_NOT_SELECTED=0, MINCUT_SET_S=1, MINCUT_SET_T=2, &
        MAXFLOW_EDMOND_KARP=10, MAXFLOW_DINIC=20

    ! Named local constants
    integer, parameter :: DEFAULT_ECAPACITY = 10, DEFAULT_VCAPACITY = 5
    integer, parameter :: MAP_NULL = -1, NOT_INITIALIZED = -1
    integer, parameter :: INTEGER_MOLD(0) = [integer ::]

    integer(I1B), parameter :: VERTEX_HANDLE_TYPE = 1_I1B, &
        EDGE_HANDLE_TYPE = 2_I1B, INVALID_HANDLE_TYPE= 0_I1B

    integer, parameter :: CONCOM_LABEL_INPROGRESS=-1, &
        CONCOM_LABEL_NOTSELECTED=0


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
      procedure :: is_initialized => graph_is_initialized
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
      procedure :: maxflow_multiple => graph_maxflow_multiple
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


    interface ! submodules

! -------------------------------
! graph_smod_centrality.f90
! -------------------------------

      module subroutine graph_betweenness(this, position_cost, position_eb, &
          position_vb, is_normalized, vselector, eselector, vmask, emask)
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
      end subroutine graph_betweenness


! -------------------
! graph_smod_flow.f90
! -------------------

      module subroutine graph_mincut(this, position_weight, mincut, labels, &
          s_list, t_list, vmask, emask, vselector, eselector)
        class(graph_t), intent(in) :: this
        integer, intent(in) :: position_weight
        real(dp), intent(out) :: mincut
        integer, intent(out), allocatable, optional :: labels(:)
        type(handle_t), intent(out), allocatable, optional :: s_list(:), t_list(:)
        logical, intent(in), optional :: vmask(:), emask(:)
        procedure(is_vertex_selected), optional :: vselector
        procedure(is_edge_selected), optional :: eselector
!
! Compute the global minimum cut of an undirected weighted graph.
!
! The routine partitions the selected vertices into two disjoint subsets
! (S,T) such that the total weight of edges crossing the partition is
! minimized. Edge weights are taken from "edges/rpar(position_weight)".
!
! Vertices and edges can be restricted using optional selector functions or
! logical array masks. If no selectors are provided, all vertices and edges
! are considered. An edge contributes to the cut only if it is selected and
! both of its endpoint vertices are also selected.
!
! The partition can be returned either as an integer label array or as two
! arrays of vertex handles corresponding to the S and T subsets.
! Unselected vertices receive label MINCUT_NOT_SELECTED and are omitted from
! both handle arrays.
!
! The implementation is based on the Stoer-Wagner global minimum cut
! algorithm and therefore requires an undirected graph.
!
! INPUT
!   this            - graph structure
!   position_weight - position in the edges(:)/rpar array containing edge
!                     weights (capacities)
!   vselector       - optional function selecting vertices participating in
!                     the computation
!   eselector       - optional function selecting edges participating in the
!                     computation
!   vmask, emask    - optional arrays selecting vertices and edges
!                     participating in the computation
!
! OUTPUT
!   mincut          - weight of the minimum cut
!   labels          - optional partition labels for all vertices:
!                       MINCUT_NOT_SELECTED - vertex not selected
!                       MINCUT_SET_S        - vertex belongs to subset S
!                       MINCUT_SET_T        - vertex belongs to subset T
!   s_list          - optional handles of vertices belonging to subset S
!   t_list          - optional handles of vertices belonging to subset T
!
      end subroutine graph_mincut


      module subroutine graph_maxflow(this, source, sink, position_capacity, &
          flow, position_mincutlabel, position_flow, &
          vmask, emask, vselector, eselector, algorithm_maxflow)
        class(graph_t), intent(inout) :: this
        type(handle_t), intent(in) :: source, sink
        integer, intent(in) :: position_capacity
        real(dp), intent(out) :: flow
        integer, intent(in), optional :: position_mincutlabel
        integer, intent(in), optional :: position_flow
        logical, intent(in), optional :: vmask(:), emask(:)
        procedure(is_vertex_selected), optional :: vselector
        procedure(is_edge_selected), optional :: eselector
        integer, intent(in), optional :: algorithm_maxflow
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
!   algorithm_maxflow - which algorithm to use:
!                         MAXFLOW_EDMOND_KARP, or
!                         MAXFLOW_DINIC
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
      end subroutine graph_maxflow


      module subroutine graph_maxflow_multiple(this, sources, sinks, &
          position_capacity, flow, position_mincutlabel, position_flow, &
          vmask, emask, vselector, eselector, algorithm_maxflow)
        class(graph_t), intent(inout) :: this
        type(handle_t), intent(in) :: sources(:), sinks(:)
        integer, intent(in) :: position_capacity
        real(dp), intent(out) :: flow
        integer, intent(in), optional :: position_mincutlabel
        integer, intent(in), optional :: position_flow
        logical, intent(in), optional :: vmask(:), emask(:)
        procedure(is_vertex_selected), optional :: vselector
        procedure(is_edge_selected), optional :: eselector
        integer, intent(in), optional :: algorithm_maxflow
!
! Maximum flow using multiple sources and sinks
!
!   algorithm_maxflow - which algorithm to use:
!                         MAXFLOW_EDMOND_KARP, or
!                         MAXFLOW_DINIC
      end subroutine graph_maxflow_multiple

    end interface ! submodules

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

      if (.not. this%is_initialized()) then
        write(fid,'("Graph not initialized (WARNING)")')
        return
      end if

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
    subroutine graph_connected_components(this, labels, lab_count, &
        position_label, vselector, eselector, vmask, emask)
      class(graph_t), intent(inout) :: this
      integer, intent(out), allocatable, optional :: labels(:)
      integer, intent(out), optional :: lab_count
      integer, intent(in), optional :: position_label
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
      logical, intent(in), optional :: vmask(:), emask(:)
!
! Identify connected components of an undirected graph or weakly connected
! components of a directed graph.
!
! The routine assigns consecutive integer labels starting from 1 to each
! connected component.
!
! Vertices and edges can be restricted using optional selector functions or
! logical array masks. If no selectors are provided, all vertices and edges are
! considered. An edge contributes to connectivity only if it is selected and
! both of its endpoint vertices are also selected.
!
! By providing or omitting optional arguments "labels" and "position_label",
! user can:
!   - obtain an array of labels (unselected vertices are marked by zero),
! or/and
!   - let the subroutine to store the labels "vertices/ipar" array. Labels of
! unselected vertices in graph are left untouched.
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
!   labels         - component label or zero for unselected vertices (optional)
!   lab_count      - optional number of connected components identified
!   this           - for selected vertices, component labels are stored to
!                    "vertices/ipar(position_label)" array
!
      integer :: i, j, iedge, idst, lab_current
      type(iterator_t) :: iterator
      type(stack_t) :: stack
      integer, allocatable :: labels0(:)
      logical, allocatable :: vmask0(:), emask0(:)
      type(adjlist_t), allocatable :: incomming_ngbs(:)

      call graph_build_selection_masks(this, vmask0, emask0, &
          vselector=vselector, eselector=eselector, vmask_provided=vmask, &
          emask_provided=emask)

      ! To identify weakly connected components in directed graphs,
      ! we need to traverse edges in the reverse direction as well.
      ! Generate adjacency lists of incoming edges that are open.
      if (this%is_directed_graph) then
        allocate(incomming_ngbs(this%nvertices))
        do i=1, this%nvertices
          call incomming_ngbs(i)%initialize()
        end do
        do iedge=1, this%nedges
          if (.not. emask0(iedge)) cycle
          idst = get_index_from_handle(this, this%edges(iedge)%dst_handle)
          if (idst==MAP_NULL) error stop 'graph_connected_component -&
              & graph contains edge with invalid destination point'
          call incomming_ngbs(idst)%add(iedge)
        end do
      end if

      ! Initialize "labels" for selected vertices
      allocate(labels0(this%nvertices))
      where (vmask0)
        labels0 = CONCOM_LABEL_INPROGRESS
      else where
        labels0 = CONCOM_LABEL_NOTSELECTED
      end where

      ! Identified components counter
      lab_current = 0

      ! Stack for deep-first graph traversal (DFS)
      call stack%initialize(chunksize=size(transfer(i,INTEGER_MOLD)))

      MAIN_LOOP: do i=1, this%nvertices
        ! Find the next unprocessed vertex and add it to the empty stack
        if (.not. vmask0(i)) cycle
        if (labels0(i) /= CONCOM_LABEL_INPROGRESS) cycle
        lab_current = lab_current + 1
        labels0(i) = lab_current
        call stack%push(transfer(i,INTEGER_MOLD))

        ! Process the stack and propagate "lab_current"
        DFS_LOOP: do while (.not. stack%empty())
          j = transfer(stack%pop(), j)

          ! Label and add allowed neighbours to the stack
          iterator = iterator_t()
          do while (this%vertices(j)%ngbs%has_next(iterator))
            call this%vertices(j)%ngbs%next(iterator, iedge)
            if (.not. emask0(iedge)) cycle
            call follow_edge(this, iedge, j, lab_current, vmask0, labels0, stack)
          end do

          ! Label and add inbound neighbours (directed graphs only)
          if (this%is_directed_graph) then
            iterator = iterator_t()
            do while (incomming_ngbs(j)%has_next(iterator))
              call incomming_ngbs(j)%next(iterator, iedge)
              call follow_edge(this, iedge, j, lab_current, vmask0, labels0, stack)
            end do
          end if

        end do DFS_LOOP

      end do MAIN_LOOP

      ! Write labels to vertex/ipar array
      if (present(position_label)) then
        where (vmask0) &
            this%vertices(1:this%nvertices)%ipar(position_label) = labels0
      end if

      if (present(labels)) call move_alloc(labels0, labels)

      if (present(lab_count)) lab_count = lab_current

      ! Clean-up (explicitly deallocating array of adjlist_t)
      if (this%is_directed_graph) deallocate(incomming_ngbs)

    end subroutine graph_connected_components


    subroutine follow_edge(this, iedge, from_vertex, lab_current, vmask, labels, stack)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: iedge, from_vertex, lab_current
      logical, intent(in) :: vmask(:)
      integer, intent(inout) :: labels(:)
      type(stack_t), intent(inout) :: stack
!
! Follow edge during labeling connected components traversal.
!
      integer :: next_vertex

      next_vertex = other_vertex_id(this, iedge, from_vertex)
      ! Assert "idst" is selected as edges to unselected vertices should
      ! have been unselected by "build_selection_masks"
      if (.not. vmask(next_vertex)) error stop &
          'graph_connected_components - selected edge has an unselected endpoint'
      ! According to its label, a neighbor "idst" can be:
      !   - unvisited selected vertex: assign current component and push to stack
      !   - already assigned to this component: ignore
      ! Other labels indicate an inconsistent graph traversal.
      if (labels(next_vertex) == CONCOM_LABEL_INPROGRESS) then
        labels(next_vertex) = lab_current
        call stack%push(transfer(next_vertex,INTEGER_MOLD))
      else if (labels(next_vertex) == lab_current) then
        continue
      else
        error stop 'graph_connected_components -&
            & neighbour belongs to another component (internal error)'
      end if
    end subroutine follow_edge


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
