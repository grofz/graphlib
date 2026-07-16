  module graph_mod
    use iso_fortran_env, only : dp => real64, i1b => int8, i8b => int64
    use conts_mod, only : queue_t, stack_t, pqueue_t, pqueue_handle_t=>handle_t, &
      PQUEUE_MIN
    use graph_adjlist_mod, only : adjlist_t, iterator_t
    implicit none (type, external)
    private

    ! Sadly, parametrized derived type (PDT) not working reliably with compilers.
    ! To avoid PDT, lets hard-code the array sizes required for the actual implementation
    integer, parameter, public :: NIV_PARS = 1, NRV_PARS = 7, NIE_PARS = 1, NRE_PARS = 2
    ! Legend
    ! V/IPAR = [type]
    ! V/RPAR = [radius, x, y, z, ?, ?, ?]
    ! E/IPAR = [type]
    ! E/RPAR = [cost, ?]

    ! Other constants
    integer, parameter :: DEFAULT_ECAPACITY = 10, DEFAULT_VCAPACITY = 5
    integer, parameter :: MAP_NULL = -1, NOT_INITIALIZED = -1
    integer, parameter :: INTEGER_MOLD(0) = [integer ::]

    integer(i1b), parameter :: &
      VERTEX_HANDLE_TYPE = 1_i1b, EDGE_HANDLE_TYPE = 2_i1b, GENERAL_HANDLE_TYPE= 0_i1b

    type, public :: handle_t
      private
      integer :: index_to_map = MAP_NULL
      integer :: version = 1
      integer(i1b) :: handle_type = GENERAL_HANDLE_TYPE
    contains
      procedure, private :: handle_eq
      generic :: operator(==) => handle_eq
    end type handle_t

    type, public :: vertex_t
      integer  :: ipar(NIV_PARS)
      real(dp) :: rpar(NRV_PARS)
      type(adjlist_t) :: ngbs ! list of outgoing edge ids
      type(handle_t) :: handle
    end type vertex_t

    type, public :: edge_t
      type(handle_t) :: src_handle, dst_handle
      integer  :: ipar(NIE_PARS)
      real(dp) :: rpar(NRE_PARS)
      type(handle_t) :: handle
    contains
      procedure :: vertex_indices => edge_vertex_indices
    end type

    type, public :: graph_t
      type(vertex_t), allocatable :: vertices(:)
      type(edge_t), allocatable :: edges(:)
      integer, allocatable :: vmap(:), emap(:)
          ! storing position of vertices/edges in "vertices"/"edges" arrays
      integer :: nvertices=NOT_INITIALIZED, nedges
      logical :: is_directed_graph=.false.
        ! .true. = edges are "one-way"
        ! .false. = edges are bi-directional
      type(queue_t) :: free_vhandles, free_ehandles
      integer :: niv, nrv, nie, nre
        ! store size of "ipar" and "rpar" arrays in vertices and edges
    contains
      procedure :: initialize => graph_initialize
      procedure :: add_vertex => graph_add_vertex
      procedure :: add_edge   => graph_add_edge
      procedure :: remove_vertex => graph_remove_vertex
      procedure :: remove_edge => graph_remove_edge
      procedure :: print => graph_print
      procedure :: connected_components => graph_connected_components
      procedure :: shortest_path => graph_shortest_path
      procedure :: maxflow => graph_maxflow
      procedure :: betweenness => graph_betweenness
    end type graph_t


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
        error stop 'borrow_handle: unknown handle_type'
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
        error stop 'return_handle: unknown handle_type'
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

      this%niv = NIV_PARS
      this%nrv = NRV_PARS
      this%nie = NIE_PARS
      this%nre = NRE_PARS
    end subroutine graph_initialize


    pure logical function graph_is_initialized(this)
      class(graph_t), intent(in) :: this
      graph_is_initialized = this%nvertices /= NOT_INITIALIZED
    end function graph_is_initialized


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
      else if (size(ipar) /= NIV_PARS .or. size(rpar) /= NRV_PARS) then
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

      integer :: ivertex

      if (handle%handle_type /= VERTEX_HANDLE_TYPE) &
          error stop 'graph_remove_vertex - invalid handle type'
      ivertex = get_index_from_handle(this, handle)
      if (ivertex == MAP_NULL) &
          error stop 'graph_remove_vertex - vertex no longer exists'

      ! All outgoing edges will be also automatically removed
!goto 111
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
!111 continue

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
      else if (size(ipar) /= NIE_PARS .or. size(rpar) /= NRE_PARS) then
        error stop 'graph_add_edge - invalid argument arrays size'
      else if (src%handle_type /= VERTEX_HANDLE_TYPE .or. dst%handle_type /= VERTEX_HANDLE_TYPE) then
        error stop 'graph_add_edge - invalid src or dst handles type'
      end if

      block
        integer :: isrc, idst

        isrc = get_index_from_handle(this, src)
        idst = get_index_from_handle(this, dst)
        if (isrc==MAP_NULL .or. idst==MAP_NULL) then
          error stop 'graph_add_edge - vertex not exists (invalid handle)'
        end if
        ! check if connection already exists
        if (find_edge_id(this, isrc, idst) /= MAP_NULL) then
          error stop 'graph_add_edge - connection already exists'
        else if (.not. this%is_directed_graph) then
          if (find_edge_id(this, idst, isrc) /= MAP_NULL) then
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

        ! add the new edge to the outgoing edges list
        call this%vertices(isrc)%ngbs%add(this%nedges)
        if (.not. this%is_directed_graph) then
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
print '("Edge ",i0,"--",i0," removed")', isrc, idst
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


    subroutine graph_print(this, fid)
      class(graph_t), intent(in) :: this
      integer, intent(in) :: fid

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
        if (this%niv>0) write(fid,*) this%vertices(i)%ipar
        if (this%nrv>0) write(fid,*) this%vertices(i)%rpar
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
        if (this%nie>0) write(fid,*) this%edges(i)%ipar
        if (this%nre>0) write(fid,*) this%edges(i)%rpar
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


    pure function find_edge_id(this, ia, ib) result(id)
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
    end function find_edge_id


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
        position_mincutlabel, position_flow, vselector, eselector)
      class(graph_t), intent(inout) :: this
      type(handle_t), intent(in) :: source, sink
      integer, intent(in) :: position_capacity
      real(dp), intent(out) :: flow
      integer, intent(in), optional :: position_mincutlabel
      integer, intent(in), optional :: position_flow
      procedure(is_vertex_selected), optional :: vselector
      procedure(is_edge_selected), optional :: eselector
!
! Maximum flow from the source to sink.
!
! INPUT
!   this              - the graph (vertex/edge data updated)
!   source            - handle to the source vertex
!   sink              - handle to the sink vertex
!   position_capacity - "edges/rpar" array item giving the edge capacity
!   vselector         - user function to select open verices (OPTIONAL)
!   eselector         - user functoin to select open edges (OPTIONAL)
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
!   position_flow        - (OPTIONAL)"edges/rpar" array item to save flow along
!                          the edge
!
      integer, parameter :: SOURCE_REACHABLE=1, SINK_REACHABLE=2, CLOSED=0
      real(dp), allocatable :: forward_capacity(:), backward_capacity(:)
      integer, allocatable :: prev_edge(:), pair_edge(:)
      integer :: source_id, sink_id
      real(dp) :: additional_flow
      type(stack_t) :: added_edges
      logical, allocatable :: vmask(:), emask(:)

      ! Set up working arrays
      allocate(forward_capacity(this%nedges), backward_capacity(this%nedges))
        ! Remaining capacity for forward and backward flow, backward_capacity
        ! is used for undirected graphs only. For directed graphs reverse
        ! edges are temporarily added to the graph.
      allocate(prev_edge(this%nvertices))
        ! Keep track to the incoming edge is.

      ! Select open edges and vertices
      call graph_build_selection_masks(this, vmask, emask, vselector=vselector, eselector=eselector)

      block
        ! Capacity of all open edges (and edges connecting open vertices) is set
        ! to their capacity given in "edges/rpar" array.
        ! Capacity of closed edges set to zero (to disabling them)
        integer :: i, ia, ib
        do i=1, this%nedges
          forward_capacity(i) = 0.0_dp
          backward_capacity(i) = 0.0_dp
          if (emask(i)) then
            ! these checks are being done in "build_selection_masks" and can be removed
            ! after testing TODO
            ia = get_index_from_handle(this, this%edges(i)%src_handle)
            ib = get_index_from_handle(this, this%edges(i)%dst_handle)
            ! orphaned edges will be silently ignored
            if (ia==MAP_NULL .or. ib==MAP_NULL) then
              error stop 'graph_maxflow - selection mask builder not working ok1'
              cycle
            end if
            if (.not. (vmask(ia) .and. vmask(ib))) then
              error stop 'graph_maxflow - selection mask builder not working ok2'
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
        else if (.not. vmask(source_id)) then
          error stop 'graph_max_flow - source is not open'
        else if (.not. vmask(sink_id)) then
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

      ! For directed graphs, reverse edges must be added to the graph. This means
      ! that "forward_capacity" will be reallocated, "backward_capacity" will no
      ! longer be needed. Reference "pair_edge" will be used instead.
      if (this%is_directed_graph) then
        block
          type(handle_t) :: edge
          integer :: nreverse_edges, i, ireverse
          real(dp), allocatable :: tmp_forward_capacity(:)

          ! count edges with non-zero capacity
          nreverse_edges = count(forward_capacity > 0.0_dp)

          ! for each non-zero capacity edge, a reverse edge is added
          allocate(tmp_forward_capacity(this%nedges+nreverse_edges), source=0.0_dp)
          allocate(pair_edge(this%nedges+nreverse_edges))
          do i=1,this%nedges
            if (.not. (forward_capacity(i)>0.0_dp)) cycle
            associate(e=>this%edges(i))
              edge = this%add_edge(e%dst_handle, e%src_handle, e%ipar, e%rpar)
            end associate
            ireverse = get_index_from_handle(this, edge)
            call added_edges%push(transfer(edge,INTEGER_MOLD))
            pair_edge(i) = ireverse
            pair_edge(ireverse) = i
            ! The capacity of reverse edges is initially stored to zero as
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

      ! The main loop of Edmonds-Karp
      ! Augment flow as long as path with non-zero capacity exists
      flow = 0.0_dp
      do
        ! Find shortest path using edges with non-zero remaining capacity
        call bfs_shortest_path(this, forward_capacity, backward_capacity, &
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

      ! Make minimum cut partition (optional)
      if (present(position_mincutlabel)) then
        block
          integer :: i
          call bfs_shortest_path(this, forward_capacity, backward_capacity, &
              source_id, sink_id, prev_edge)
          do i=1,this%nvertices
            associate(label=>this%vertices(i)%ipar(position_mincutlabel))
              if (.not. vmask(i)) then
                label = CLOSED
              else if (i==source_id) then
                label = SOURCE_REACHABLE
              else if (prev_edge(i)/=MAP_NULL) then
                label = SOURCE_REACHABLE
              else
                label = SINK_REACHABLE
                ! NOTE - actually we do not know if this vertex had any connection
                !        with sink. Additional check for isolated vertices may be
                !        added later (TODO)
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


    subroutine bfs_shortest_path(this, forward_capacity, backward_capacity, &
        source_id, sink_id, prev_edge)
      class(graph_t), intent(in) :: this
      real(dp), intent(in) :: forward_capacity(:), backward_capacity(:)
      integer, intent(in) :: source_id, sink_id
      integer, intent(out) :: prev_edge(:)
!
! BFS from source to sink using only edges with non-zero capacity.
! The path (if exists) can be tracked using "prev_edge" array.
! Undirected graphs: traversing edge as SRC-<DST uses forward_capacity,
! traversing edge as DST->SRC uses backward_capacity.
! Directed graphs: backward_capacity is not used.
!
      type(queue_t) :: q
      integer :: current_id, iedge, ngb_id
      type(iterator_t) :: iterator

      call q%initialize(chunksize=size(transfer(current_id,INTEGER_MOLD)))
      call q%enqueue(transfer(source_id,INTEGER_MOLD))
      prev_edge = MAP_NULL

      do while(.not. q%empty() .and. prev_edge(sink_id)==MAP_NULL)
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
      ! Now it is possible use "prev_edge(sink_id)" to see if path from
      ! source to sink was found and back-track the path back to source.
    end subroutine bfs_shortest_path


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
! Back-track the path and
!  (i) find the bottleneck remaining capacity, or
!  (ii) update remaining capacity along the path.
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
      if (current_id /= source_id) error stop 'process_path - could not reach source'
    end subroutine process_path


! TODO make a version for unweighted graps
! - position_cost to be optional
! - using "dist_dijkstra" real array, and "dist_bfs" integer arrray
! - factor out the djiksta code to special procudeure
! - write out special procedure for BFS search
    ! -------------------------------
    ! Betweenness (Brandes algorithm)
    ! -------------------------------
    subroutine graph_betweenness(this, position_cost, position_eb, position_vb, &
        is_normalized, vselector, eselector, vmask, emask)
      class(graph_t), intent(inout) :: this
      integer, intent(in) :: position_cost
      integer, intent(in) :: position_eb, position_vb
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
      integer :: i, ia, ib, id_s
      real(dp), allocatable :: dist(:)
        ! shortest distance from "s" to "v"
      type(stack_t), allocatable :: prev(:)
        ! list of immediate predecessors of "v" on shortest paths
      integer(I8B), allocatable :: sigma(:)
        ! the number of unique shortest paths from "s" to "v"
      type(pqueue_t) :: pqueue
      type(stack_t) :: stack
      type(pqueue_handle_t), allocatable :: phas(:)
        ! handles to vertices added to Dijksta's priority queue
      logical, allocatable :: visited(:)
        ! denote visited nodes during Dijksta search
      real(dp), allocatable :: delta(:)
        ! dependency of the source on vertex "v"
      integer, parameter :: DIST_SHORTER=1, DIST_SAME=0, DIST_LONGER=-1


      ! Mark selected edges and count selected_vertices
      call graph_build_selection_masks(this, vmask0, emask0, &
          vselector=vselector, eselector=eselector, &
          vmask_provided=vmask, emask_provided=emask)

      ! Local Dijkstra arrays
      allocate(dist(this%nvertices), prev(this%nvertices), sigma(this%nvertices))
      allocate(phas(this%nvertices), visited(this%nvertices))
      do i=1, this%nvertices
        call prev(i)%initialize(size(transfer(id_s,INTEGER_MOLD)))
      end do
      call pqueue%initialize(chunksize=size(transfer(id_s,INTEGER_MOLD)))
      call stack%initialize(chunksize=size(transfer(id_s,INTEGER_MOLD)))
      allocate(delta(this%nvertices))

      ! Initialize edge and vertex betweeness
      where (emask0) &
          this%edges(1:this%nedges)%rpar(position_eb) = 0.0_dp
      where (vmask0) &
          this%vertices(1:this%nvertices)%rpar(position_vb) = 0.0_dp

      ! Main loop over all source vectors
      SRC_LOOP: do id_s=1, this%nvertices
        if (.not. vmask0(id_s)) cycle
if (mod(id_s,500)==0) print '("Source is ",i0," out of ",i0)', id_s, this%nvertices

        ! Initialize Dijkstra's search from source
        dist = huge(dist)
        do i=1, this%nvertices
          call prev(i)%clear()
          phas(i) = pqueue_handle_t()
        end do
        sigma = 0_I8B
        visited = .false.

        dist(id_s) = 0.0_dp
        sigma(id_s) = 1_I8B
        phas(id_s) = pqueue%insert(transfer(id_s,INTEGER_MOLD), dist(id_s))

        ! STEP 1 - Modified Dijkstra loop
        DJIKSTRA_LOOP: do while (.not. pqueue%empty())
          block
            integer :: id_v, id_u, edge_uv
            real(dp) :: dist_to_v
            type(iterator_t) :: iterator
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
          end block
        end do DJIKSTRA_LOOP

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
              associate(c=>this%edges(iedge)%rpar(position_eb))
                c = c + delta_edge
              end associate
              ! pass it to a predecessor node
              delta(id_u) = delta(id_u) + delta_edge
            end do
            if (id_v /= id_s) then
              associate(c=>this%vertices(id_v)%rpar(position_vb))
                c = c + delta(id_v)
              end associate
            end if
          end block
        end do BACKPASS_LOOP

      end do SRC_LOOP

      ! Divide by two for undirected graphs (all paths were counted twice)
      if (.not. this%is_directed_graph) then
        associate(eb=>this%edges(1:this%nedges)%rpar(position_eb))
          where (emask0) eb = 0.5_dp * eb
        end associate
        associate(vb=>this%vertices(1:this%nvertices)%rpar(position_vb))
          where (vmask0) vb = 0.5_dp * vb
        end associate
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
          if (n>=3) then
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
          if (n>=2) then
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
      deallocate(phas) ! got run time error (compiler bug?), explicit deallocation solved this

    contains
      integer function compare_dist(new, old)
        real(dp), intent(in) :: old, new
        real(dp), parameter :: REL_TOL = 1.0e5 * epsilon(1.0_dp)
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
    end subroutine graph_betweenness


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

  end module graph_mod
