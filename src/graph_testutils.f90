module graph_testutils_mod
  use graph_mod, only : graph_t, handle_t
  use graph_user_mod
  use iso_fortran_env, only : dp=>real64, output_unit
  use parse_mod, only : string_t, split_nonempty
  implicit none (type, external)

  real(dp), parameter :: INVALID = -11.0_dp

  type testsample_t
    type(graph_t) :: g
    real(dp) :: expected_mincut = INVALID
    integer, allocatable :: expected_s(:), expected_t(:)
  end type

contains

  subroutine parse_lines(lines, current_line, ts)
    type(string_t), intent(in) :: lines(:)
    integer, intent(inout) :: current_line
    type(testsample_t), intent(out) :: ts
!
! TODO Unpolished, hasted implementation
!
    type(string_t), allocatable :: words(:)
    integer :: ios, i

    do
      ! split current line into tokens
      call split_nonempty(lines(current_line)%str, ' ', words)
      current_line = current_line + 1
      if (size(words)==0) exit ! current_line is empty line 
      if (words(1)%str(1:1)=='#') then
        print '(a)', '<'//lines(current_line-1)%str//'>'
        cycle ! current_line is a comment
      end if

      ! identify first token
      select case(words(1)%str)
      case ('EXPECTED_MINCUT')
        read(words(2)%str,*,iostat=ios) ts%expected_mincut
        if (ios/=0) then
          print '(a)', lines(current_line-1)%str
          error stop 'error reading number'
        end if
      case ('EXPECTED_SET_A')
        if (allocated(ts%expected_s)) deallocate(ts%expected_s)
        allocate(ts%expected_s(size(words)-1))
        do i=2, size(words)
          read(words(i)%str,*,iostat=ios) ts%expected_s(i-1)
          if (ios/=0) then
            print '(a)', lines(current_line-1)%str
            print '(a)', '<'//words(i)%str//'>'
            error stop 'error reading numbers'
          end if
        end do
      case ('EXPECTED_SET_B')
        if (allocated(ts%expected_t)) deallocate(ts%expected_t)
        allocate(ts%expected_t(size(words)-1))
        do i=2, size(words)
          read(words(i)%str,*,iostat=ios) ts%expected_t(i-1)
          if (ios/=0) then
            print '(a)', lines(current_line-1)%str
            print '(a)', '<'//words(i)%str//'>'
            error stop 'error reading numbers'
          end if
        end do
      case ('EDGES')
        block
          integer, allocatable :: ia(:), ib(:), v_ipars(:,:), e_ipars(:,:), cons(:,:)
          real(dp), allocatable :: weights(:), v_rpars(:,:), e_rpars(:,:)
          integer :: k, i1, i2, ios1, ios2, nv, ne
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
            read(tokens(3)%str,*,iostat=ios2) w1
            if (ios1/=0 .or. ios2/=0) then
              print '(a)', lines(current_line+k)%str
              print '(a)', '<'//words(i)%str//'>'
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
          nv = maxval(ia)
          nv = max(nv, maxval(ib))
          ne = size(ia)
          allocate(v_ipars(VSIZE_IPAR, nv), v_rpars(VSIZE_RPAR, nv))
          allocate(e_ipars(ESIZE_IPAR, ne), e_rpars(ESIZE_RPAR, ne))
          e_rpars(EPOS_WEIGHT,:) = weights(:)
          allocate(cons(2,ne))
          cons(1,:) = ia
          cons(2,:) = ib
          call graph_from_arrays(ts%g, cons, erdata=e_rpars)
        end block

      case default
        print '(a)', '<'//words(1)%str//'> is unknown token'
        error stop
      end select

      if (current_line > size(lines)) exit ! end of page reached
    end do

   !call ts%g%print(output_unit)

  end subroutine parse_lines


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
    block
      character(len=:), allocatable :: str
      str = 'undirected'
      if (g%is_directed()) str = 'directed'
      print '("new ",a," graph with ",i0," vertices and ",i0," edges")', &
          str, g%nvertices, g%nedges
    end block
  end subroutine graph_from_arrays
end module graph_testutils_mod