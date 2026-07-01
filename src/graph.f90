  module graph_mod
    use iso_fortran_env, only : dp => real64, i1b => int8
    use conts_mod, only : queue_t
    implicit none (type, external)
    private

    ! Parametrized derived types not working with gfortran
    ! - to avoid PDT, we must hard-code the array sizes required for the actual implementation
    integer, parameter, public :: NIV_PARS = 0, NRV_PARS = 0, NIE_PARS = 0, NRE_PARS = 0

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
      integer, allocatable :: ngbs(:) ! list of outgoing edges id
      type(handle_t) :: handle
    end type

    type edge_t
      type(handle_t) :: src_handle, dst_handle
      integer  :: ipar(NIE_PARS)
      real(dp) :: rpar(NRE_PARS)
      type(handle_t) :: handle
    end type

    type, public :: graph_t
      type(vertex_t), allocatable :: vertices(:)
      type(edge_t), allocatable :: edges(:)
      integer, allocatable :: vmap(:), emap(:)
          ! to store the position of vertices/edges in "vertices"/"edges" arrays
      integer :: nvertices=NOT_INITIALIZED, nedges
      logical :: is_directed_graph=.false.
        ! .true. = edges are "one-way"
        ! .false. = edges are bi-directional
      type(queue_t) :: free_vhandles, free_ehandles
    end type

  contains

      pure function handle_eq(a, b) result(eq)
        class(handle_t), intent(in) :: a, b
        logical :: eq
        eq = a%version==b%version .and. &
             a%index_to_map==b%index_to_map .and. &
             a%handle_type==b%handle_type
      end function handle_eq


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
        integer :: id
!
! Return position of a vertex/edge in array using handle. If handle refers to
! the vertex/edge no longer in array, MAP_NULL value is returned
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


      function graph_add_vertex(this, ipar, rpar) result(handle)
        class(graph_t), intent(inout) :: this
        integer, intent(in) :: ipar(:)
        real, intent(in) :: rpar(:)
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
          allocate (new_vertex%ngbs(0))
          new_vertex%handle = handle
        end associate
        this%vmap(handle%index_to_map) = this%nvertices
      end function graph_add_vertex


      subroutine reallocate_vertex(this, handle, newid)
        class(graph_t), intent(inout) :: this
        type(handle_t), intent(in) :: handle
        integer, intent(in) :: newid

        integer :: oldid

        if (handle%handle_type /= VERTEX_HANDLE_TYPE) &
            error stop 'reallocate_vertex - wrong handle type' 
        oldid = get_index_from_handle(this, handle)
        if (oldid == MAP_NULL) &
            error stop 'reallocate_vertex - vertex no more exists'
        if (newid < 1 .or. newid > this%nvertices) &
            error stop 'reallocate_vertex - newid out of bounds'

        ! copy vertex and update record in "vmap"
        this%vertices(newid) = this%vertices(oldid)
        this%vmap(handle%index_to_map) = newid
      end subroutine reallocate_vertex


      function graph_add_edge(this, src, dst, ipar, rpar) result(handle)
        class(graph_t), intent(inout) :: this
        type(handle_t), intent(in) :: src, dst
        integer, intent(in) :: ipar(:)
        integer, intent(in) :: rpar(:)
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

          ! include new edge to the outgoing edges list
          ! TODO more effective, without reallocation?
          this%vertices(isrc)%ngbs = [this%vertices(isrc)%ngbs, this%nedges]
          if (.not. this%is_directed_graph) then
            this%vertices(idst)%ngbs = [this%vertices(idst)%ngbs, this%nedges]
          end if
        end block
      end function graph_add_edge


      ! TODO
      ! subroutine reallocate_edge


      function list_of_ngbs(this, isrc) result(idsts)
        class(graph_t), intent(in) :: this
        integer, intent(in) :: isrc
        integer :: idsts(size(this%vertices(isrc)%ngbs))
!
! Return an array of neighbors of "isrc" vertex.
! For directed graph, only outbound neighbours are listed.
!
        integer :: i, ia, ib

        associate (ngbs => this%vertices(isrc)%ngbs)
          do i=1, size(ngbs)
            if (ngbs(i) <= 0 .or. ngbs(i) > this%nedges) then
              error stop 'list_of_ngbs - item in ngbs is out of bounds'
            end if
            ia = get_index_from_handle(this, this%edges(ngbs(i))%src_handle)
            ib = get_index_from_handle(this, this%edges(ngbs(i))%dst_handle)
            if (ia==MAP_NULL .or. ib==MAP_NULL) then
              error stop 'list_of_ngbs - edge has invalid handles (vertex no more exists)'
            else if (ia/=isrc .and. ib/=isrc) then
              error stop 'list_of_ngbs - no edge endpoint is "isrc" (should not happen)'
            end if

            if (.not. this%is_directed_graph) then ! bi-directional graph
              if (ia==isrc) then
                idsts(i) = ib
              else if (ib==isrc) then
                idsts(i) = ia
              else
                error stop 'list_of_ngbs - should not be reachable'
              end if
            else ! directed graph
              if (ia==isrc) then
                idsts(i) = ib
              else
                error stop 'list_of_ngbs - wrong source in directed graph'
              end if
            end if
          end do
        end associate
      end function list_of_ngbs


      function get_connection_index(this, ia, ib) result(id)
        class(graph_t), intent(in) :: this
        integer, intent(in) :: ia, ib
        integer :: id

        integer :: list(size(this%vertices(ia)%ngbs))
        integer :: i

        id = MAP_NULL
        list = list_of_ngbs(this, ia)
        do i=1, size(list)
          if (list(i) /= ib) cycle
          id = this%vertices(ia)%ngbs(i)
          exit
        end do
      end function get_connection_index

  end module graph_mod
