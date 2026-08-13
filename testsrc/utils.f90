module testutils_mod
  use iso_fortran_env, only : dp=>real64, output_unit
  use graph_mod, only : graph_t, handle_t
  use map_mod
  use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR
  use parse_mod, only : string_t, split_nonempty
  implicit none (type, external)
  private
  public parse_lines, graph_from_arrays

  real(dp), parameter :: INVALID = -11.0_dp

  type, public :: testsample_t
    type(graph_t) :: g
    real(dp) :: expected_mincut = INVALID
    real(dp) :: expected_maxflow(2) = INVALID
    integer, allocatable :: expected_s(:), expected_t(:), expected_scc(:)
    real(dp), allocatable :: positions(:,:)
    integer, allocatable :: sources(:), sinks(:)
    real(dp), allocatable :: expected_eb(:), expected_vb(:)
    logical :: is_directed_graph = .false.
  contains
    procedure :: reset => testsample_reset
  end type

contains

  subroutine testsample_reset(this)
    class(testsample_t), intent(inout) :: this

    this%expected_mincut = INVALID
    this%expected_maxflow = INVALID
    if (allocated(this%expected_s)) deallocate(this%expected_s)
    if (allocated(this%expected_t)) deallocate(this%expected_t)
    if (allocated(this%expected_scc)) deallocate(this%expected_scc)
    if (allocated(this%positions)) deallocate(this%positions)
    if (allocated(this%sources)) deallocate(this%sources)
    if (allocated(this%sinks)) deallocate(this%sinks)
    if (allocated(this%expected_eb)) deallocate(this%expected_eb)
    if (allocated(this%expected_vb)) deallocate(this%expected_vb)
!TODO keep it commented for now, test_maxflow fails otherwise
!   this % is_directed_graph = .false.
    call this%g%initialize()
  end subroutine testsample_reset


  subroutine parse_lines(lines, current_line, ts)
    type(string_t), intent(in) :: lines(:)
    integer, intent(inout) :: current_line
    type(testsample_t), intent(inout) :: ts
!
! IN
!   lines        - strings to be parsed
!
! IN/OUT
!   current_line - points to the next line to be processed
!   ts           - parse data to the structure
!
! TODO Unpolished, hasted implementation
!
    type(string_t), allocatable :: words(:)

    call ts%reset()
    do
      ! split current line into tokens
      current_line = current_line + 1
      call split_nonempty(lines(current_line-1)%str, ' ', words)
      if (size(words)==0 .and. ts%g%nvertices>0) then
        ! empty line between different graphs
        exit
      else if (size(words)==0) then
        ! empty line separating different sections
        cycle
      else if (words(1)%str(1:1)=='#') then
        ! a comments
#ifdef DEBUG
        print '(a)', '<'//lines(current_line-1)%str//'>'
#endif
        cycle
      end if

      ! identify first token of current_line
      select case(words(1)%str)
      case ('EXPECTED_MINCUT')
        block
          real(dp), allocatable :: arr(:)
          call parse_list_of_reals(words, arr)
          if (size(arr)==1) then
            ts%expected_mincut = arr(1)
          else
            print '(a)', lines(current_line-1)%str
            error stop 'error reading expected mincut'
          end if
        end block

      case ('EXPECTED_MAXFLOW')
        block
          real(dp), allocatable :: arr(:)
          call parse_list_of_reals(words, arr)
          if (size(arr)==2) then
            ts%expected_maxflow(1:2) = arr
          else
            print '(a)', lines(current_line-1)%str
            error stop 'error reading expected maxflow (two values expected)'
          end if
        end block

      case('EXPECTED_VB')
        call parse_list_of_reals(words, ts%expected_vb)

      case('EXPECTED_EB')
        call parse_list_of_reals(words, ts%expected_eb)

      case ('EXPECTED_SET_A')
        call parse_list_of_ints(words, ts%expected_s)

      case ('EXPECTED_SET_B')
        call parse_list_of_ints(words, ts%expected_t)

      case ('EXPECTED_SCC')
        call parse_list_of_ints(words, ts%expected_scc)

      case('SOURCES')
        call parse_list_of_ints(words, ts%sources)

      case('SINKS')
        call parse_list_of_ints(words, ts%sinks)

      case('DIRECTED')
        ts%is_directed_graph = .true.

      case('VERTICES')
        ! must be before 'EDGES'
        call parse_vertices(lines, current_line, ts)

      case ('EDGES')
        call parse_edges(lines, current_line, ts)

      case default
        print '(a)', '<'//words(1)%str//'> is unknown token'
        error stop
      end select

      if (current_line > size(lines)) exit ! end of file reached
    end do

  end subroutine parse_lines


  subroutine parse_vertices(lines, current_line, ts)
    type(string_t), intent(in) :: lines(:)
    integer, intent(inout) :: current_line
    type(testsample_t), intent(inout) :: ts

    real(dp), allocatable :: xyz(:)
    integer :: k, ios(3)
    real(dp) :: x, y, z
    type(string_t), allocatable :: tokens(:)
    allocate(xyz(0))
    k = 0
    do
      call split_nonempty(lines(current_line+k)%str, ' ', tokens)
      if (size(tokens)==0) then
        exit
      else if (size(tokens)/=3) then
        print '(a)', lines(current_line+k)%str
        error stop 'parse_vertices - three postions expected'
      end if
      read(tokens(1)%str,*,iostat=ios(1)) x
      read(tokens(2)%str,*,iostat=ios(2)) y
      read(tokens(3)%str,*,iostat=ios(3)) z
      if (any(ios/=0)) then
        print '(a)', lines(current_line+k)%str
        error stop 'error parse vertices'
      end if
      xyz = [xyz, x, y, z]

      k = k+1
      if (current_line+k>size(lines)) exit
    end do
    current_line = current_line+k

    ! save "position" component
    ts%positions=reshape(xyz, shape=[3, size(xyz)/3])
  end subroutine parse_vertices


  subroutine parse_edges(lines, current_line, ts)
    type(string_t), intent(in) :: lines(:)
    integer, intent(inout) :: current_line
    type(testsample_t), intent(inout) :: ts

    integer, allocatable :: ia(:), ib(:), v_ipars(:,:), e_ipars(:,:), cons(:,:)
    real(dp), allocatable :: weights(:), v_rpars(:,:), e_rpars(:,:)
    integer :: k, i1, i2, ios1, ios2, nv, ne, ios
    real(dp) :: w1
    type(string_t), allocatable :: tokens(:)
    allocate(ia(0),ib(0),weights(0))
    k = 0
    do
      call split_nonempty(lines(current_line+k)%str, ' ', tokens)
      if (size(tokens)==0) then
        exit
      end if
      read(tokens(1)%str,*,iostat=ios) i1
      if (ios /= 0) exit
      read(tokens(2)%str,*,iostat=ios1) i2
      if (size(tokens)>=3) then
        read(tokens(3)%str,*,iostat=ios2) w1
      else
        ios2 = 0
        w1 = 1.0_dp
      end if
      if (ios1/=0 .or. ios2/=0) then
        print '(a)', lines(current_line+k)%str
        error stop 'error reading numbers'
      end if
      ia = [ia, i1]
      ib = [ib, i2]
      weights = [weights, w1]

      k = k+1
      if (current_line+k>size(lines)) exit
    end do
    current_line = current_line+k

    ! prepare graph arrays
    if (allocated(ts%positions)) then
      nv = size(ts%positions,2)
    else
      nv = maxval(ia)
      nv = max(nv, maxval(ib))
    end if
    ne = size(ia)
    allocate(v_ipars(VSIZE_IPAR, nv), v_rpars(VSIZE_RPAR, nv))
    allocate(e_ipars(ESIZE_IPAR, ne), e_rpars(ESIZE_RPAR, ne))
    e_rpars(EPOS_WEIGHT,:) = weights(:)
    if (allocated(ts%positions))then
      v_rpars(VPOS_X+0,:) = ts%positions(1,:)
      v_rpars(VPOS_X+1,:) = ts%positions(2,:)
      v_rpars(VPOS_X+2,:) = ts%positions(3,:)
    end if
    allocate(cons(2,ne))
    cons(1,:) = ia
    cons(2,:) = ib
    call graph_from_arrays(ts%g, cons, erdata=e_rpars, vrdata=v_rpars, &
        is_directed_graph=ts%is_directed_graph)
  end subroutine parse_edges


  subroutine parse_list_of_ints(words, numbers)
    type(string_t), intent(in) :: words(:)
    integer, allocatable, intent(out) :: numbers(:)

    integer :: i, ios

    if (allocated(numbers)) deallocate(numbers)
    allocate(numbers(size(words)-1))
    do i=2, size(words)
      read(words(i)%str,*,iostat=ios) numbers(i-1)
      if (ios/=0) then
        print '(a," is not a number ")', '<'//words(i)%str//'>'
        error stop 'parse_list_of_ints fail'
      end if
    end do
  end subroutine parse_list_of_ints


  subroutine parse_list_of_reals(words, numbers)
    type(string_t), intent(in) :: words(:)
    real(dp), allocatable, intent(out) :: numbers(:)

    integer :: i, ios

    if (allocated(numbers)) deallocate(numbers)
    allocate(numbers(size(words)-1))
    do i=2, size(words)
      read(words(i)%str,*,iostat=ios) numbers(i-1)
      if (ios/=0) then
        print '(a," is not a number ")', '<'//words(i)%str//'>'
        error stop 'parse_list_of_reals fail'
      end if
    end do
  end subroutine parse_list_of_reals


  subroutine graph_from_arrays(g, cons, vidata, vrdata, eidata, erdata, &
      is_directed_graph)
    type(graph_t), intent(inout) :: g
    integer, intent(in) :: cons(:,:)
    integer, intent(in), optional :: vidata(:,:), eidata(:,:)
    real(dp), intent(in), optional :: vrdata(:,:), erdata(:,:)
    logical, intent(in), optional :: is_directed_graph
!
! TODO
!
    integer :: nv, ne, i
    integer :: v_idat(VSIZE_IPAR), e_idat(ESIZE_IPAR)
    real(dp) :: v_rdat(VSIZE_RPAR), e_rdat(ESIZE_RPAR)
    type(handle_t), allocatable :: vhandles(:)
    type(handle_t) :: handle

    ! initialize by arbitrary values
    v_idat = 0
    e_idat = 0
    v_rdat = 0.0_dp
    e_rdat = 0.0_dp

    ! verify correct array dimensions
    nv = maxval(cons)
    ne = size(cons, dim=2)
    if (size(cons,1)/=2) error stop &
        'graph_from_arrays - first dimension of cons must be two'
    if (present(vidata)) then
      if (size(vidata,1)/=VSIZE_IPAR .or. size(vidata,2)/=nv) error stop &
        'graph_from_arrays - vidata has invalid shape'
    end if
    if (present(vrdata)) then
      if (size(vrdata,1)/=VSIZE_RPAR .or. size(vrdata,2)/=nv) error stop &
        'graph_from_arrays - vrdata has invalid shape'
    end if
    if (present(eidata)) then
      if (size(eidata,1)/=ESIZE_IPAR .or. size(eidata,2)/=ne) error stop &
        'graph_from_arrays - eidata has invalid shape'
    end if
    if (present(erdata)) then
      if (size(erdata,1)/=ESIZE_RPAR .or. size(erdata,2)/=ne) error stop &
        'graph_from_arrays - erdata has invalid shape'
    end if
    if (any(cons<1)) error stop &
        'graph_from_arrays - all integers in cons must be positive'

    ! Initialize the graph
    call g%initialize(is_directed_graph=is_directed_graph)

    ! Add vertices
    allocate(vhandles(nv))
    do i=1, nv
      if (present(vidata)) v_idat = vidata(:,i)
      if (present(vrdata)) v_rdat = vrdata(:,i)
      vhandles(i) = g%add_vertex(v_idat, v_rdat)
    end do

    ! Add edges
    do i=1, ne
      if (present(eidata)) e_idat = eidata(:,i)
      if (present(erdata)) e_rdat = erdata(:,i)
      handle = g%add_edge(vhandles(cons(1,i)), vhandles(cons(2,i)), e_idat, e_rdat)
    end do

    ! report
#ifdef DEBUG
    block
      character(len=:), allocatable :: str
      str = 'undirected'
      if (g%is_directed()) str = 'directed'
      print '("new ",a," graph with ",i0," vertices and ",i0," edges")', &
          str, g%nvertices, g%nedges
    end block
#endif
  end subroutine graph_from_arrays
end module testutils_mod
