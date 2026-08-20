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


  submodule (graph_mod) graph_smod_connectivity
!
! Label connected components.
! Tarjan's algorithm (determine strongly connected components).
! Khanns algorithm (determine topological levels)
!
    implicit none (type, external)

  contains

! ----------------------------------------------------------------------------
!   subroutine graph_connected_components(this, labels, lab_count, &
!       position_label, vselector, eselector, vmask, emask)
! ----------------------------------------------------------------------------
    module procedure graph_connected_components
      integer :: i, j, iedge, idst, lab_current
      type(iterator_t) :: iterator
      type(stack_t) :: stack
      integer, allocatable :: labels0(:)
      logical, allocatable :: vmask0(:), emask0(:)
      type(adjlist_t), allocatable :: incomming_ngbs(:)

      ! Verify position_label flag is not out of bounds
      if (present(position_label)) then
        if (position_label < 1 .or. position_label > VSIZE_IPAR) error stop &
            'graph_connected_components - position_label out of bounds'
      end if

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
          idst = this%index_from_handle(this%edges(iedge)%dst_handle)
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
      call stack%initialize(chunksize=size(transfer(1,INTEGER_MOLD)))

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

    end procedure graph_connected_components


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


! ----------------------------------------------------------------------------
!   subroutine graph_strongly_connected_components(this, labels, lab_count, &
!       position_label, vselector, eselector, vmask, emask)
! ----------------------------------------------------------------------------
    module procedure graph_strongly_connected_components
      logical, allocatable :: vmask0(:), emask0(:)
      integer :: scc_counter, id_counter
      integer, allocatable :: components(:), discovered_id(:), lowlink(:)
      type(iterator_t), allocatable :: iterators(:)
      type(stack_t) :: dfs_stack, scc_stack
      integer :: istart, iu, ingb, iedge
      integer, parameter :: UNASSIGNED=0, UNVISITED=0

      if (.not. this%is_directed_graph) error stop &
          'graph_scc: algorithm requires a directed graph'

      ! Verify position_label flag is not out of bounds
      if (present(position_label)) then
        if (position_label < 1 .or. position_label > VSIZE_IPAR) error stop &
            'graph_scc - position_label out of bounds'
      end if

      call graph_build_selection_masks(this, vmask0, emask0, &
          vselector=vselector, eselector=eselector, vmask_provided=vmask, &
          emask_provided=emask)

      ! Initialize working arrays and global counters
      scc_counter = 0
      id_counter = 0
      allocate(components(this%nvertices), source=UNASSIGNED)
        ! vertices not yet assigned to a particular SCC
      allocate(discovered_id(this%nvertices), source=UNVISITED)
        ! vertices not yet discovered
      allocate(lowlink(this%nvertices))
        ! lowest DFS discovery index reachable from this vertex
      allocate(iterators(this%nvertices), source=iterator_t())
        ! to keep track of explored outgoing edges

      ! Initialize stacks to store vertice positions in graph array
      call dfs_stack%initialize(chunksize=size(transfer(1,INTEGER_MOLD)))
      call scc_stack%initialize(chunksize=size(transfer(1,INTEGER_MOLD)))

      ! Loop to make sure all vertices are discovered
      MAIN_LOOP: do istart=1, this%nvertices
        ! Push an open, undiscovered vertex as a first item on the stack
        if ((.not. vmask0(istart)) .or. discovered_id(istart)/=UNVISITED) cycle
        call dfs_stack%push(transfer(istart,INTEGER_MOLD))

        DFS_LOOP: do while(.not. dfs_stack%empty())
          iu = transfer(dfs_stack%peek(), iu)
          ! If U discovered, initialize it and add it on SCC stack
          if (discovered_id(iu)==UNVISITED) then
            id_counter = id_counter+1
            discovered_id(iu) = id_counter
            lowlink(iu) = discovered_id(iu)
            call scc_stack%push(transfer(iu,INTEGER_MOLD))
          end if

          ! Check if U has remaining unexplored neighbours
          if (this%vertices(iu)%ngbs%has_next(iterators(iu))) then
            call this%vertices(iu)%ngbs%next(iterators(iu), iedge)
            if (.not. emask0(iedge)) cycle DFS_LOOP
            ingb = other_vertex_id(this, iedge, iu)

            if (discovered_id(ingb)==UNVISITED) then
              ! Case A: unvisited neighbour (down one level of recursion)
              call dfs_stack%push(transfer(ingb,INTEGER_MOLD))
            else if (components(ingb)==UNASSIGNED) then
              ! Case B: neighbour currently on the SCC stack
              ! Update lowlink using its DFS index
              lowlink(iu) = min(lowlink(iu), discovered_id(ingb))
            else
              ! Case C: neighbour already assigned to an SCC
              continue ! ignore this neighbour
            end if
            cycle DFS_LOOP
          end if

          ! All neighbours of U have been processed
          block
            integer :: itmp
            itmp = transfer(dfs_stack%pop(),itmp)
              ! up one level of recursion
          end block

          ! DFS returns to the parent. Propagate the minimum reachable
          ! DFS index upwards.
          if (.not. dfs_stack%empty()) then
            block
              integer :: iparent
              iparent = transfer(dfs_stack%peek(),iparent)
              lowlink(iparent) = min(lowlink(iparent), lowlink(iu))
            end block
          end if

          ! Check if U is an SCC root
          if (lowlink(iu)==discovered_id(iu)) then
            ! New SCC is complete...
            scc_counter = scc_counter + 1
            ! ...label and remove all SCC vertices from SCC stack
            block
              integer :: imoved
              do
                imoved = transfer(scc_stack%pop(),imoved)
                components(imoved) = scc_counter
                if (imoved==iu) exit
              end do
            end block
          end if

        end do DFS_LOOP

      end do MAIN_LOOP

      if (.not. scc_stack%empty()) error stop &
          'graph_scc - unprocessed vertices in stack (internal error)'

      ! Write labels to vertex/ipar array
      if (present(position_label)) then
        where (vmask0) &
            this%vertices(1:this%nvertices)%ipar(position_label) = components
      end if

      if (present(labels)) call move_alloc(components, labels)

      if (present(lab_count)) lab_count = scc_counter

      ! Clean-up (explicitly deallocate array of iterator_t)
      deallocate(iterators)

    end procedure graph_strongly_connected_components


! ----------------------------------------------------------------------------
!   subroutine graph_topological_levels(this, levels, components, vmask, &
!       emask, vselector, eselector)
! ----------------------------------------------------------------------------
    module procedure graph_topological_levels
      logical, allocatable :: vmask0(:), emask0(:)
      integer, allocatable :: indegree_c(:), levels_c(:), components0(:)
      type(adjlist_t), allocatable :: list_c(:)
      integer :: icomp, vertex_id, iedge, ngb_id
      integer, parameter :: UNASSIGNED=0, UNDEFINED=-1
      type(queue_t) :: queue
      type(iterator_t) :: viterator, eiterator

      if (.not. this%is_directed_graph) error stop &
          'graph_topological_levels: algorithm requires a directed graph'

      ! Select sub-graph
      call graph_build_selection_masks(this, vmask0, emask0, &
          vmask_provided=vmask, emask_provided=emask, vselector=vselector, &
          eselector=eselector)

      ! Build or verify components (i.e. group of vertices)
      block
        integer :: max_component_label, k

        if (present(components)) then
          ! Assert valid components provided
          if (size(components)/=this%nvertices) error stop &
              'graph_topological_levels - size of components invalid'
          if (any(components<0)) error stop &
              'graph_topological_levels - negative value in components'
          if ((any(components/=UNASSIGNED .and. .not. vmask0)) .or. &
              (any(components==UNASSIGNED .and. vmask0))) error stop &
              'graph_topological_levels - selected vertices list mismatch with assigned components list'
          components0 = components
        else
          ! Make each vertex its own component.
          allocate(components0(this%nvertices), source=UNASSIGNED)
          k = 1
          do vertex_id=1, this%nvertices
            if (.not. vmask0(vertex_id)) cycle
            components0(vertex_id) = k
            k = k + 1
          end do
        end if

        max_component_label = maxval(components0)
        allocate(indegree_c(max_component_label), source=0)
        allocate(levels_c(max_component_label), source=UNDEFINED)
        allocate(list_c(max_component_label))
        do icomp=1, max_component_label
          call list_c(icomp)%initialize()
        end do
        do vertex_id=1, this%nvertices
          if (.not. vmask0(vertex_id)) cycle
          call list_c(components0(vertex_id))%add(vertex_id)
        end do
      end block

      ! Calculate indegree, i.e. the number of incomming edges to each component
      ! (group of vertices)
      block
        integer :: src, dst
        do iedge=1, this%nedges
          if (.not. emask0(iedge)) cycle
          src = this%index_from_handle(this%edges(iedge)%src_handle)
          dst = this%index_from_handle(this%edges(iedge)%dst_handle)
          ! ignore internal edges within a component
          if (components0(src)==components0(dst)) cycle
          indegree_c(components0(dst)) = indegree_c(components0(dst))+1
        end do
      end block

      ! Initialize resulting array
      allocate(levels(this%nvertices), source=0)
        ! unselected vertices will keep this initial value (0)

      ! Enqueue the components with zero indegree
      call queue%initialize(chunksize=size(transfer(1,INTEGER_MOLD)))
      do icomp=1, size(indegree_c)
        if (indegree_c(icomp)>0) cycle
        levels_c(icomp) = 1 ! starting topological level
        if (list_c(icomp)%size()/=0) then
          call queue%enqueue(transfer(icomp,INTEGER_MOLD))
        else
          ! If items in provided "components" array do not form a serie of
          ! consecutive values, groups of zero vertices may be present.
          ! To make them harmless, We assign the initial level to these groups
          ! but will not add then to the queue.
          continue
        end if
      end do

      ! Process components in topological order
      QLOOP: do while (.not. queue%empty())
        icomp = transfer(queue%dequeue(),icomp)
        viterator = iterator_t()
        ! Loop over every vertex in the "icomp" group...
        VLOOP: do while (list_c(icomp)%has_next(viterator))
          call list_c(icomp)%next(viterator, vertex_id)
          ! ... assign level to this vertex
          levels(vertex_id) = levels_c(icomp)
          ! ... and loop over its neighbours and decrease indegree
          !     of the destination group
          eiterator = iterator_t()
          NGBLOOP: do while (this%vertices(vertex_id)%ngbs%has_next(eiterator))
            call this%vertices(vertex_id)%ngbs%next(eiterator, iedge)
            if (.not. emask0(iedge)) cycle
            ngb_id = other_vertex_id(this, iedge, vertex_id)
            ! ignore intra-group vertices
            if (components0(ngb_id)==components0(vertex_id)) cycle

            associate (X=>indegree_c(components0(ngb_id)))
              if (X<1) error stop 'graph_topological_levels -&
                  & indegree drops under zero (internal error)'
              X = X - 1
              if (X==0) then
                ! No more dependencies for the group, add it to the queue
                if (levels_c(components0(ngb_id))/=UNDEFINED) error stop &
                    & 'graph_topological_levels - rediscovering group (internal error)'
                levels_c(components0(ngb_id)) = levels_c(icomp)+1
                call queue%enqueue(transfer(components0(ngb_id),INTEGER_MOLD))
              end if
            end associate
          end do NGBLOOP
        end do VLOOP
      end do QLOOP

      ! Verify all groups were processed.
      ! Unprocessed groups mean that a cycle is present in the graph.
      if (any(levels_c==UNDEFINED)) error stop 'graph_topological_levels -&
          & graph contains directed cycles, topological levels cannot be derermined'

    end procedure graph_topological_levels


! ----------------------------------------------------------------------------
!   module function graph_verify_topological_levels(this, levels, &
!       components, vmask, emask, vselector, eselector) result(is_ok)
! ----------------------------------------------------------------------------
    module procedure graph_verify_topological_levels
      logical, allocatable :: vmask0(:), emask0(:)
      integer :: iedge, iu, iv

      if (.not. this%is_directed_graph) error stop &
          'graph_verify_topological_levels - a directed graph required'

      if (size(levels)/=this%nvertices) error stop &
          'graph_verify_topological_levels - levels size is invalid'

      if (present(components)) then
        if (size(components)/=this%nvertices) error stop &
          'graph_verify_topological_levels - components size is invalid'
      end if

      call graph_build_selection_masks(this, vmask0, emask0, &
          vmask_provided=vmask, emask_provided=emask, &
          vselector=vselector, eselector=eselector)

      is_ok = .false. ! must reach end if all is ok

      ! level>0 is expected for selected vertices
      ! level=0 is expected for unselected vertices
      if (any(levels<0)) return
      if (any(levels==0 .and. vmask0)) return
      if (any(levels>0 .and. .not. vmask0)) return

      do iedge=1, this%nedges
        if (.not. emask0(iedge)) cycle
        iu = this%index_from_handle(this%edges(iedge)%src_handle)
        iv = this%index_from_handle(this%edges(iedge)%dst_handle)
        if (present(components)) then
          if (components(iv)==components(iu)) then
            ! intra-component edge
            ! all vertices within a component must have same level
            if (levels(iu)/=levels(iv)) return
            cycle
          end if
        end if
        ! verify topological level definition is met
        if (levels(iu)>=levels(iv)) return
      end do
      is_ok = .true.
    end procedure graph_verify_topological_levels

  end submodule graph_smod_connectivity
