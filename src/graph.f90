  module graph_mod
    use iso_fortran_env, only : dp => real64, i1b => int8
    use conts_mod, only : queue_t
    use graph_adjlist_mod, only : adjlist_t, iterator_t
    implicit none (type, external)
    private

    ! Sadly, parametrized derived type (PDT) not working reliably with compilers.
    ! To avoid PDT, lets hard-code the array sizes required for the actual implementation
    integer, parameter, public :: NIV_PARS = 1, NRV_PARS = 4, NIE_PARS = 1, NRE_PARS = 0

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

    type vertex_t
      integer  :: ipar(NIV_PARS)
      real(dp) :: rpar(NRV_PARS)
      type(adjlist_t) :: ngbs ! list of outgoing edges ids
      type(handle_t) :: handle
    end type

    type edge_t
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
    end type

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
        if (this%free_vhandles%size()==0) call graph_increase_vertices_capacity(this)
        if (this%free_vhandles%size()==0) error stop 'borrow_handle - no more handles available V'
        handle = transfer(this%free_vhandles%dequeue(), handle)
      case(EDGE_HANDLE_TYPE)
        if (this%free_ehandles%size()==0) call graph_increase_edges_capacity(this)
        if (this%free_ehandles%size()==0) error stop 'borrow_handle - no more handles available E'
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
! the vertex/edge that is no longer in array, MAP_NULL is returned
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
        call graph_increase_vertices_capacity(this, new_capacity)
        new_capacity = DEFAULT_ECAPACITY
        if (present(ecapacity)) new_capacity = ecapacity
        call graph_increase_edges_capacity(this, new_capacity)
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


    subroutine graph_increase_vertices_capacity(this, new_capacity)
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
    end subroutine graph_increase_vertices_capacity


    subroutine graph_increase_edges_capacity(this, new_capacity)
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
    end subroutine graph_increase_edges_capacity


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
      block
        type(iterator_t) :: iterator
        integer :: iedge
        iterator = iterator_t()
        do while(this%vertices(ivertex)%ngbs%has_next(iterator))
          call this%vertices(ivertex)%ngbs%next(iterator, iedge)
          call graph_remove_edge(this, this%edges(iedge)%handle)
          iterator = iterator_t()
        end do
      end block

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
        if (get_connection_index(this, isrc, idst) /= MAP_NULL) then
          error stop 'graph_add_edge - connection already exists'
        else if (.not. this%is_directed_graph) then
          if (get_connection_index(this, idst, isrc) /= MAP_NULL) then
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


      function list_of_ngbs(this, isrc) result(idsts)
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

          idst = edge_other_vertex_id(this, iedge, isrc)
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


      function edge_other_vertex_id(this, iedge, ia) result(ib)
        class(graph_t), intent(in) :: this
        integer, intent(in) :: iedge, ia
        integer ib

        integer :: i1, i2

        i1 = get_index_from_handle(this, this%edges(iedge)%src_handle)
        i2 = get_index_from_handle(this, this%edges(iedge)%dst_handle)
        if (i1==MAP_NULL .or. i2==MAP_NULL) then
          error stop 'edge_other_vertex_id - edge has invalid handles (vertex no more exists)'
        else if (ia==i1) then
          ib = i2
        else if (ia==i2) then
          ib = i1
        else
          error stop 'edge_other_vertex_id - no edge endpoint is "ia" (should not happen)'
        end if
      end function edge_other_vertex_id


      function get_connection_index(this, ia, ib) result(id)
        class(graph_t), intent(in) :: this
        integer, intent(in) :: ia, ib
        integer :: id

        integer :: iedge
        type(iterator_t) :: iterator

        id = MAP_NULL
        iterator = iterator_t()
        do while (this%vertices(ia)%ngbs%has_next(iterator))
          call this%vertices(ia)%ngbs%next(iterator, iedge)
          if (edge_other_vertex_id(this, iedge, ia) == ib) then
            id = iedge
            exit
          end if
        end do
      end function get_connection_index


      function edge_vertex_indices(this, graph) result(ids)
        class(edge_t), intent(in) :: this
        type(graph_t), intent(in) :: graph
        integer ::ids(2)

        ids(1) = get_index_from_handle(graph, this%src_handle)
        ids(2) = get_index_from_handle(graph, this%dst_handle)
      end function edge_vertex_indices

  end module graph_mod
