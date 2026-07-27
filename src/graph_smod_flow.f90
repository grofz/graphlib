  submodule(graph_mod) graph_smod_flow
!
! Stoer-Wagner algorithm for min-cut partition
! Edmond-Karp and Dinics algorithms for maximum flow
!
    implicit none (type, external)

  contains

    ! ============
    ! Stoer-Wagner
    ! ============

! -----------------------------------------------------------------------------
!   subroutine graph_mincut(this, position_weight, mincut, labels, &
!       s_list, t_list,  vmask, emask, vselector, eselector)
! -----------------------------------------------------------------------------
    module procedure graph_mincut
      type(graph_t) :: g
      type(handle_t) :: handle
      type(stack_t), allocatable :: contracted_vertices(:)
      type(stack_t) :: mincut_vertices
      integer, allocatable :: labels0(:)
      type(pqueue_handle_t), allocatable :: pqueue_handles(:)
      type(pqueue_t) :: scoreboard

      ! Stoer-Wagner algorithm works for undirected graphs only
      if (this%is_directed_graph) error stop &
          'graph_mincut - algorithm requires an undirected graph only'

      ! For s_list and t_list optional arguments, either both or none must
      ! be provided
      if (present(s_list) .neqv. present(t_list)) error stop &
          'graph_mincut - s_list/t_list arguments: both or none must be provided'

      ! Prepare the working graph and the array of stacks with handles to the
      ! vertices in the original graph.
      block
        integer :: i, k
        type(handle_t), allocatable :: vertices0(:)

        call this%copy(g, vselector=vselector, eselector=eselector, &
            vmask=vmask, emask=emask, new_vertices=vertices0)
        if (g%nvertices < 2) error stop &
            'graph_mincut - at least two vertices are required'
        allocate(contracted_vertices(g%nvertices))
        do i=1, g%nvertices
          if (g%vertices(i)%handle%index_to_map /= i) error stop &
            'graph_mincut - handle index does not match array index (internal error)'
          call contracted_vertices(i)%initialize( &
              chunksize=size(transfer(handle,INTEGER_MOLD)))
        end do
        allocate(labels0(this%nvertices), source=MINCUT_NOT_SELECTED)
        i = 1
        do k=1, size(vertices0)
          if (get_index_from_handle(g, vertices0(k))==MAP_NULL) cycle
          ! The k-th vertex in the original graph was selected and now is the
          ! i-th vertex in the working graph. Push the handle from the original
          ! graph as the first item to the stack in the working graph
          call contracted_vertices(i)%push( &
              transfer(this%vertices(k)%handle,INTEGER_MOLD))
          ! Label vertices in original graph LAB_SET_S or LAB_NOT_SELECTED.
          ! LAB_SET_S label is temporary, some vertices will be relabeled when
          ! the min-cut is determined
          labels0(k) = MINCUT_SET_S
          i = i+1
        end do
        if (i-1/=g%nvertices) error stop 'graph_mincut -&
            & valid items count mismatch working graph size (internal error)'
      end block

      ! The priority queue stores the handles of vertices not yet moved
      ! to the set A, the vertex connectivity with vertices already present
      ! in the set A is stored as the priority in the queue.
      call scoreboard%initialize( &
          chunksize=size(transfer(handle,INTEGER_MOLD)), ordering=PQUEUE_MAX)
      allocate(pqueue_handles(g%nvertices))

      ! Initialize global mincut weight and corresponding vertices list
      call mincut_vertices%initialize(chunksize=size(transfer(handle,INTEGER_MOLD)))
      mincut = huge(mincut)

      block
        type(handle_t) :: s_handle, t_handle
        real(dp) :: cut_weight
        MAIN_LOOP: do
          ! Phase 1
          call find_st(g, position_weight, s_handle, t_handle, cut_weight, &
              scoreboard, pqueue_handles)

          ! Update min-cut weight and vertices list if a lower weight found
          if (cut_weight < mincut) then
            mincut = cut_weight
            mincut_vertices = contracted_vertices(t_handle%index_to_map)
          end if

          ! Phase 2 - contract vertices S and T
          if (g%nvertices <= 2) exit MAIN_LOOP
          call contract_st(g, position_weight, s_handle, t_handle, contracted_vertices)
        end do MAIN_LOOP
      end block

      ! Consume mincut_vertices to label vertices from the winning mincut
      do while (.not. mincut_vertices%empty())
        handle = transfer(mincut_vertices%pop(), handle)
        associate(lab=>labels0(get_index_from_handle(this, handle)))
          if (lab /= MINCUT_SET_S) error stop &
              'graph_mincut - label0 value unexpected (internal error)'
          lab = MINCUT_SET_T
        end associate
      end do

      ! Make s_list and t_list (if required by user)
      if (present(s_list) .and. present(t_list)) then
        block
          integer :: is, it, k
          allocate(s_list(count(labels0==MINCUT_SET_S)))
          allocate(t_list(count(labels0==MINCUT_SET_T)))
          is = 1
          it = 1
          do k=1, size(labels0)
            if (labels0(k)==MINCUT_SET_S) then
              s_list(is) = this%vertices(k)%handle
              is = is+1
            else if (labels0(k)==MINCUT_SET_T) then
              t_list(it) = this%vertices(k)%handle
              it = it+1
            end if
          end do
          if (is-1/=size(s_list) .or. it-1/=size(t_list)) error stop &
              'graph_mincut - miscalculation while building s/t_lists (internal error)'
        end block
      end if

      ! Return labels if required by user
      if (present(labels)) call move_alloc(labels0, labels)

      ! Clean-up (explicitly deallocate array of pqueue_handle_t)
      deallocate(pqueue_handles)

    end procedure graph_mincut


    subroutine find_st(g, position_weight, s_handle, t_handle, cut_weight, &
          scoreboard, phandles)
      type(graph_t), intent(in) :: g
      integer, intent(in) :: position_weight
      type(handle_t), intent(out) :: s_handle, t_handle
      real(dp), intent(out) :: cut_weight
      type(pqueue_t), intent(inout) :: scoreboard
      type(pqueue_handle_t), intent(inout) :: phandles(:)
!
! Perform one Stoer-Wagner minimum-cut phase.
!
! Vertices are iteratively added to the growing set A. On return,
! s_handle and t_handle are the last two vertices added to A, and
! cut_weight is the weight of the corresponding s-t cut.
!
! The working priority queue contains vertices assumed currently outside
! of set A.
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
      end do

      ! the queue contains the last vertex T and its connectivity with
      ! all remaining vertices
      t_handle = transfer(scoreboard%pop(top_priority=cut_weight), t_handle)

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


    ! ===================
    ! Edmond-Karp / Dinic
    ! ===================

! -----------------------------------------------------------------------------
!   subroutine graph_maxflow(this, source, sink, position_capacity, &
!       flow, position_mincutlabel, position_flow, &
!       vmask, emask, vselector, eselector, algorithm_maxflow)
! -----------------------------------------------------------------------------
    module procedure graph_maxflow
      integer, parameter :: &
        CLOSED=0, SOURCE_REACHABLE=1, SINK_REACHABLE=2, DISCONNECTED=3
      integer, parameter :: NOT_DISCONNECTED=-1
      real(dp), allocatable :: forward_capacity(:), backward_capacity(:)
      integer, allocatable :: prev_edge(:), pair_edge(:)
      integer :: source_id, sink_id
      type(stack_t) :: added_edges
      logical, allocatable :: vmask0(:), emask0(:)
      integer :: algorithm_maxflow0

      ! Select algorithm - Edmond-Karp or Dinic
      algorithm_maxflow0 = MAXFLOW_DINIC ! default value
      if (present(algorithm_maxflow)) algorithm_maxflow0 = algorithm_maxflow

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

      ! The core of the algorithm
      select case(algorithm_maxflow0)
      case(MAXFLOW_EDMOND_KARP)
        call edmond_karp_loop(this, forward_capacity, backward_capacity, &
            source_id, sink_id, prev_edge, pair_edge, flow, position_flow)
      case(MAXFLOW_DINIC)
        call dinic_loop(this, forward_capacity, backward_capacity, &
            source_id, sink_id, prev_edge, pair_edge, flow, position_flow)
      case default
        error stop 'graph_maxflow - invalid algorihm id'
      end select

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

    end procedure graph_maxflow


    subroutine edmond_karp_loop(this, forward_capacity, backward_capacity, &
        source_id, sink_id, prev_edge, pair_edge, flow, position_flow)
      class(graph_t), intent(inout) :: this
      real(dp), intent(inout) :: forward_capacity(:), backward_capacity(:)
      integer, intent(in) :: source_id, sink_id, pair_edge(:)
      integer, intent(inout) :: prev_edge(:)
      real(dp), intent(out) :: flow
      integer, intent(in), optional :: position_flow
!
! The main loop of Edmonds-Karp
! Augment flow as long as path with non-zero capacity exists
!
      real(dp) :: additional_flow

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
    end subroutine edmond_karp_loop


    subroutine dinic_loop(this, forward_capacity, backward_capacity, &
        source_id, sink_id, prev_edge, pair_edge, flow, position_flow)
      class(graph_t), intent(inout) :: this
      real(dp), intent(inout) :: forward_capacity(:), backward_capacity(:)
      integer, intent(in) :: source_id, sink_id, pair_edge(:)
      integer, intent(inout) :: prev_edge(:)
      real(dp), intent(out) :: flow
      integer, intent(in), optional :: position_flow
!
! The main loop of Dinic's algorithm
! Augment flow as long as path with non-zero capacity exists
!
      real(dp) :: additional_flow
      integer, allocatable :: levels(:)
      integer :: i
      type(iterator_t), allocatable :: iterators(:)

      ! How many edges from the source vertex in the residual graph
      allocate(levels(this%nvertices))

      ! To save what branches were already explored
      allocate(iterators(this%nvertices))

      flow = 0.0_dp
      MAIN: do
        ! Update level by a complete BFS traversal from the source
        call bfs_residual_search(this, forward_capacity, backward_capacity, &
            source_id, 0, prev_edge, levels)

        ! No more path from source to sink with non-zero capacity exists
        if (levels(sink_id)==MAP_NULL) exit MAIN
print '("Dinic: Sink distance is ",i0," levels")', levels(sink_id)

        ! Initialize iterators for DFS
        do i = 1, this%nvertices
          iterators(i) = iterator_t()
        end do

        DFSLOOP: do
          ! Find shortest path using edges with non-zero remaining capacity
          call dfs_dinic(this, levels, forward_capacity, backward_capacity, &
              source_id, sink_id, prev_edge, iterators)

          ! No more path exists with the current levels
          if (prev_edge(sink_id)==MAP_NULL) exit DFSLOOP

          ! The flow can be augmented. How much flow can we send?
          call process_path(this, forward_capacity, backward_capacity, &
              source_id, sink_id, prev_edge, pair_edge, additional_flow, .false.)
print *, 'Dinic: Current flow is ', flow,'. Augmenting by ',additional_flow,'.'
          ! Update capacity of the network
          call process_path(this, forward_capacity, backward_capacity, &
              source_id, sink_id, prev_edge, pair_edge, additional_flow, .true., &
              position_flow)
          flow = flow + additional_flow
        end do DFSLOOP

      end do MAIN

      ! Explicitly deallocate array of iterator_t
      deallocate(iterators)

    end subroutine dinic_loop


    subroutine bfs_residual_search(this, forward_capacity, backward_capacity, &
        source_id, target_id, prev_edge, levels)
      class(graph_t), intent(in) :: this
      real(dp), intent(in) :: forward_capacity(:), backward_capacity(:)
      integer, intent(in) :: source_id, target_id
      integer, intent(out) :: prev_edge(:)
      integer, intent(out), optional :: levels(:)
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
!   levels           - (optional) distance (number of edges) from source:
!                              0 - source vertex
!                            > 0 - all rachable vertices - positive integer
!                       MAP_NULL - all not reachable vertices
!
      type(queue_t) :: q
      integer :: current_id, iedge, ngb_id
      type(iterator_t) :: iterator

      call q%initialize(chunksize=size(transfer(current_id,INTEGER_MOLD)))
      call q%enqueue(transfer(source_id,INTEGER_MOLD))
      prev_edge = MAP_NULL
      if (present(levels)) then
        levels = MAP_NULL
        levels(source_id) = 0
      end if

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

          ! Save the distance from source
          if (present(levels)) then
            levels(ngb_id) = levels(current_id) + 1
          end if
        end do NGBS_LOOP
      end do
      ! Now it is possible use "prev_edge(target_id)" to see if path from
      ! source to target exists and back-track the path back to source.
      ! The source vertex remains MAP_NULL as it has no predecessor.
    end subroutine bfs_residual_search


    subroutine dfs_dinic(this, levels, forward_capacity, backward_capacity, &
              source_id, sink_id, prev_edge, iterators)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: levels(:), source_id, sink_id
      real(dp), intent(in) :: forward_capacity(:), backward_capacity(:)
      integer, intent(out) :: prev_edge(:)
      type(iterator_t), intent(inout) :: iterators(:)
!
! DFS for Dinic's algorithm
!
      type(stack_t) :: s
      integer :: current_id, iedge, ngb_id
      real(dp) :: capacity
      integer :: i

      ! Initialize prev_edge and an empty stack. Push source vertex to stack
      prev_edge = MAP_NULL
      call s%initialize(chunksize=size(transfer(current_id,INTEGER_MOLD)))
      call s%push(transfer(source_id,INTEGER_MOLD))

      STACK_LOOP: do while (.not. s%empty())
        current_id = transfer(s%peek(), current_id)

        NGBS_LOOP: do while (this%vertices(current_id)%ngbs%has_next( &
            iterators(current_id)))
          call this%vertices(current_id)%ngbs%next_noadvance( &
              iterators(current_id), iedge)
          ngb_id = other_vertex_id(this, iedge, current_id)

          ! Ignore neighbours on the same or lower BFS level. Any traversable
          ! residual edge must connect consecutive levels.
          !
          ! TODO
          ! In production version, this condition may become
          ! "if (levels(ngb_id) == levels(current_id)+1) then"
          ! to make it a bit faster and avoid following assertion
          if (levels(ngb_id) > levels(current_id)) then
            ! Determine available capacity of the "current--neighbour" edge
            if (ngb_id == get_index_from_handle(this,this%edges(iedge)%dst_handle)) then
              ! forward edge
              capacity = forward_capacity(iedge)
            else if (ngb_id == get_index_from_handle(this,this%edges(iedge)%src_handle)) then
              ! backward edge
              ! no backward edge can appear in directed graph
              if (this%is_directed_graph) error stop &
                  'dfs_dinic - assertion for directed graph fails'
              capacity = backward_capacity(iedge)
            else
              error stop 'dfs_dinic - should not reach this branch'
            end if

            if (capacity > 0.0_dp) then
              ! Internal check that levels do not differ by more than one
              ! as this would imply invalid BFS level identification.
              if (.not. levels(ngb_id)==levels(current_id)+1) &
                error stop 'dinic_dfs - following edge increses lebele by more than one (internal error)'
              ! If the edge to the neighbour has free capacity, mark which
              ! edge was used to come-in to neighbour, and then...
              prev_edge(ngb_id) = iedge
              if (ngb_id==sink_id) then
                ! ... exit if target was reached, or...
                exit STACK_LOOP
              else
                ! ... push neighbour on stack, to be explored next.
                call s%push(transfer(ngb_id,INTEGER_MOLD))
                cycle STACK_LOOP
              end if
            end if
          end if

          ! It is not possible to follow path to neighbour, try the next one.
          call this%vertices(current_id)%ngbs%advance(iterators(current_id))

        end do NGBS_LOOP

        ! All neighbours of the current vertex were explored without reaching
        ! targer. Pop current vertex from the stack, and advance iterator
        ! of the previous vertex (to not explore current vertex again)
        current_id = transfer(s%pop(), current_id)
        if (.not. s%empty()) then
          associate (p=> transfer(s%peek(), current_id))
            associate (ngbs_of_p=>this%vertices(p)%ngbs)
              if (ngbs_of_p%has_next(iterators(p))) &
                  call ngbs_of_p%advance(iterators(p))
            end associate
          end associate
        end if

      end do STACK_LOOP

      ! Now it is possible use "prev_edge(target_id)" to see if path from
      ! source to target exists and back-track the path back to source.
      ! The source vertex remains MAP_NULL as it has no predecessor.
    end subroutine dfs_dinic


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


! -----------------------------------------------------------------------------
!   subroutine graph_maxflow_multiple(this, sources, sinks, &
!       position_capacity, flow, position_mincutlabel, position_flow, &
!       vmask, emask, vselector, eselector, algorithm_maxflow)
! -----------------------------------------------------------------------------
    module procedure graph_maxflow_multiple
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
          position_flow=position_flow, vmask=vmask0, emask=emask0, &
          algorithm_maxflow=algorithm_maxflow)

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

    end procedure graph_maxflow_multiple

  end submodule graph_smod_flow
