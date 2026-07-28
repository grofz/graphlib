! ----------------------------------------------
! A unit test for grpah_betweenness calculation
! ----------------------------------------------
program betweeness
  use graph_testutils_mod, only : testsample_t, parse_lines
  use parse_mod, only : string_t, read_strings
  use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR, &
      VTUIO_MASK, &
      POS_COST=>EPOS_WEIGHT, EPOS_EB, VPOS_VB, VPOS_TYPE
  use graph_mod, only : graph_t, handle_t
  use iso_fortran_env, only : dp=>real64, output_unit
  use vtuio_mod, only : vtuio_read, vtuio_write, vtuio_data_t
  implicit none (type, external)

  type(testsample_t) :: ts
  type(string_t), allocatable :: lines(:)
  integer :: i, k
  character(len=2) :: numstr
  real(dp) :: time_start, time_end, time(2), mu, var
  real(dp), allocatable :: expected_vb(:), expected_eb(:)
  logical, allocatable :: vmask(:), emask(:)
  type(graph_t) :: g

  abstract interface
    subroutine test_graph_ai(graph, is_directed, evb, eeb)
      use graph_mod, only : graph_t, handle_t
      use iso_fortran_env, only : dp=>real64, output_unit
      implicit none (type, external)
      type(graph_t), intent(inout) :: graph
      logical, intent(in) :: is_directed
      real(dp), allocatable, intent(out) :: evb(:), eeb(:)
    end subroutine
  end interface
  interface
    pure elemental function almost_equal(a,b) result(res)
      use iso_fortran_env, only : dp=>real64
      implicit none
      real(dp), intent(in) :: a, b
      logical :: res
    end function
    subroutine statistic(x, mask, mu, var)
      import DP
      implicit none
      real(dp), intent(in) :: x(:)
      logical, intent(in) :: mask(:)
      real(dp), intent(out) :: mu, var
    end subroutine
  end interface
  procedure(test_graph_ai) :: test_graph1, test_graph2, test_graph3, test_graph4
  type(vtuio_data_t) :: vtudata

  ! betweenness tests
  lines = read_strings('assets/betweenness_sample_graphs.txt')
  i = 1
  k = 1
  call vtudata%add_item('edge_betweenness', EPOS_EB, 3, 1, 4)
  call vtudata%add_item('vertex_betweenness', VPOS_VB, 2, 1, 4)

  do while (i<= size(lines))
    ts%is_directed_graph = .false.
    call parse_lines(lines, i, ts)
    ! unweighted graph - values of vertices/edge betweenness
    call ts%g%betweenness( &
        position_eb=EPOS_EB, &
        position_vb=VPOS_VB )
    write(numstr,'(i2.2)') k
    call vtuio_write('betweennesssample'//numstr, ts%g, VTUIO_MASK, vtudata=vtudata)
    k = k+1

   !print '("Vertex central betweenness: expected/got")'
   !print '(*(g0,1x))', ts%expected_vb
   !print '(*(g0,1x))', ts%g%vertices(1:ts%g%nvertices)%rpar(VPOS_VB)
    print '("Passed? ",l1)', all(almost_equal(ts%g%vertices(1:ts%g%nvertices)%rpar(VPOS_VB), ts%expected_vb))
   !print '("Edge central betweenness: expected/got")'
   !print '(*(g0,1x))', ts%expected_eb
   !print '(*(g0,1x))', ts%g%edges(1:ts%g%nedges)%rpar(EPOS_EB)
    print '("Passed? ",l1)', all(almost_equal(ts%g%edges(1:ts%g%nedges)%rpar(EPOS_EB), ts%expected_eb))
    print *
  end do

end program


subroutine statistic(x, mask, mu, var)
  use iso_fortran_env, only : dp=>real64
  real(dp), intent(in) :: x(:)
  logical, intent(in) :: mask(:)
  real(dp), intent(out) :: mu, var

  real(dp) :: m1, m2, n

  m1 = sum(x, mask=mask)
  m2 = sum(x*x, mask=mask)
  n = real(count(mask), kind=DP)

  if (n>00) then
    mu = m1/n
    var = m2/n - mu**2
  else
    mu = -1.0_dp
    var = -1.0_dp
  end if
end subroutine


pure elemental function almost_equal(a,b) result(res)
  use iso_fortran_env, only : dp=>real64
  implicit none
  real(dp), intent(in) :: a, b
  logical :: res
  real(dp), parameter :: eps = 5.0_dp*epsilon(1.0_dp)

  res = abs(a-b) < eps * max(1.0_dp, abs(a), abs(b))
end function
