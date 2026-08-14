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


  submodule(graph_mod) graph_smod_centrality
!
! Betweenness centrality (Brandes algorithm)
!
    implicit none (type, external)

  contains

! -----------------------------------------------------------------------------
!   subroutine graph_betweenness(this, position_cost, position_eb, position_vb, &
!       is_normalized, vselector, eselector, vmask, emask)
! -----------------------------------------------------------------------------
    module procedure graph_betweenness
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

      ! Verify position flags are not out of range
      if (present(position_cost)) then
        if (position_cost < 1 .or. position_cost > ESIZE_RPAR) error stop &
            'graph_betweenness - position_cost is out of bounds'
      end if
      if (present(position_vb)) then
        if (position_vb < 1 .or. position_vb > VSIZE_RPAR) error stop &
            'graph_betweenness - position_vb is out of bounds'
      end if
      if (present(position_eb)) then
        if (position_eb < 1 .or. position_eb > ESIZE_RPAR) error stop &
            'graph_betweenness - position_eb is out of bounds'
      end if

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

    end procedure graph_betweenness


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

  end submodule graph_smod_centrality
