  submodule(graph_mod) graph_smod_flow
!
! Stoer-Wagner algorithm for min-cut partition
!
    implicit none (type, external)

  contains

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
          if (get_index_from_handle(this, vertices0(k))==MAP_NULL) cycle
          ! The k-th vertex in the original graph was selected and now is the
          ! i-th vertex in the working graph. Push the handle from the original
          ! graph as the first item to the stack in the working graph
          call contracted_vertices(i)%push(transfer(vertices0(k),INTEGER_MOLD))
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

  end submodule graph_smod_flow