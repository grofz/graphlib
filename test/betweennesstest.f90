! ----------------------------------------------
! A unit test for grpah_betweenness calculation
! ----------------------------------------------
program betweeness
  use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
  use iso_fortran_env, only : dp=>real64, output_unit
  implicit none (type, external)

  type(graph_t) :: g
  real(dp), allocatable :: expected_vb(:), expected_eb(:)
  integer :: i

  abstract interface
    subroutine test_graph_ai(g, is_directed, evb, eeb)
      use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
      use iso_fortran_env, only : dp=>real64, output_unit
      implicit none (type, external)
      type(graph_t), intent(inout) :: g
      logical, intent(in) :: is_directed
      real(dp), allocatable, intent(out) :: evb(:), eeb(:)
    end subroutine
  end interface
  procedure(test_graph_ai) :: test_graph1, test_graph2, test_graph3

  integer, parameter :: pos_cost=1, pos_eb=2, pos_vb=1

  do i=1,3
    ! Sample graphs
    select case(i)
    case(1)
      call test_graph1(g, .false., expected_vb, expected_eb)
    case(2)
      call test_graph2(g, .false., expected_vb, expected_eb)
    case(3)
      call test_graph3(g, .true., expected_vb, expected_eb)
    end select

    call g%betweenness(pos_cost, position_eb=pos_eb, position_vb=pos_vb, is_normalized=.true.)
    print '("Vertex central betweenness: expected/got")'
    print '(*(g0.6,1x))', expected_vb
    print '(*(g0.6,1x))', g%vertices(1:g%nvertices)%rpar(pos_vb)
    print '("Edge central betweenness: expected/got")'
    print '(*(g0.6,1x))', expected_eb
    print '(*(g0.6,1x))', g%edges(1:g%nedges)%rpar(pos_eb)
    print '("Passed? ",l1)', all(g%vertices(1:g%nvertices)%rpar(pos_vb)==expected_vb)
    print '("Passed? ",l1)', all(g%edges(1:g%nedges)%rpar(pos_eb)==expected_eb)
    print *
  end do

 !call g%print(output_unit)
end program


subroutine test_graph1(g, is_directed, expected_vb, expected_eb)
  use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
  use iso_fortran_env, only : dp=>real64, output_unit
  implicit none (type, external)
  type(graph_t), intent(inout) :: g
  real(dp), allocatable, intent(out) :: expected_vb(:), expected_eb(:)
  logical, intent(in) :: is_directed

  type(handle_t), allocatable :: edges(:), vertices(:)
  real(dp) :: v_rpar(NRV_PARS), e_rpar(NRE_PARS)
  integer :: v_ipar(NIV_PARS), e_ipar(NIE_PARS)

  print '("Graph 1  directed =",l1)', is_directed
  call g%initialize(is_directed_graph=is_directed)

  allocate(vertices(0:4))
  v_ipar = 1
  v_rpar = 0.0_dp
  v_rpar(4) = 0.5_dp

  v_rpar(1:3) = real([0.0, 0.0, 0.0], dp)
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([0.0, 1.0, 0.0], dp)
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([1.0, 1.0, 0.0], dp)
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([1.0, 0.0, 0.0], dp)
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([2.0, 2.0, 0.0], dp)
  vertices(4) = g%add_vertex(v_ipar, v_rpar)

  e_ipar = 10
  e_rpar = 0.0_dp
  e_rpar(1) = 1.0_dp
  allocate(edges(5))
  edges(1) = g%add_edge(vertices(0),vertices(1),e_ipar,e_rpar)
  edges(2) = g%add_edge(vertices(1),vertices(2),e_ipar,e_rpar)
  edges(3) = g%add_edge(vertices(2),vertices(3),e_ipar,e_rpar)
  edges(4) = g%add_edge(vertices(0),vertices(3),e_ipar,e_rpar)
  edges(5) = g%add_edge(vertices(2),vertices(4),e_ipar,e_rpar)

  ! undirected
  expected_vb = real([0.5,1.0,3.5,1.0,0.0],dp) / 6.0_dp
  expected_eb = real([2.5,3.5,3.5,2.5,4.0],dp) / 10.0_dp
end subroutine test_graph1


subroutine test_graph2(g, is_directed, expected_vb, expected_eb)
  use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
  use iso_fortran_env, only : dp=>real64, output_unit
  implicit none (type, external)
  type(graph_t), intent(inout) :: g
  real(dp), allocatable, intent(out) :: expected_vb(:), expected_eb(:)
  logical, intent(in) :: is_directed


  type(handle_t), allocatable :: edges(:), vertices(:)
  real(dp) :: v_rpar(NRV_PARS), e_rpar(NRE_PARS)
  integer :: v_ipar(NIV_PARS), e_ipar(NIE_PARS)

  print '("Graph 2  directed =",l1)', is_directed
  call g%initialize(is_directed_graph=is_directed)

  allocate(vertices(0:5))
  v_ipar = 1
  v_rpar = 0.0_dp
  v_rpar(4) = 0.5_dp

  v_rpar(1:3) = real([0.0, 1.0, 0.0], dp)
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([1.0, 1.0, 0.0], dp)
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([0.0, 0.0, 0.0], dp)
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([1.0, 0.0, 0.0], dp)
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([2.0, 0.0, 0.0], dp)
  vertices(4) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([2.0, 1.0, 0.0], dp)
  vertices(5) = g%add_vertex(v_ipar, v_rpar)

  e_ipar = 10
  e_rpar = 0.0_dp
  e_rpar(1) = 1.0_dp
  allocate(edges(7))
  edges(1) = g%add_edge(vertices(0),vertices(1),e_ipar,e_rpar)
  edges(2) = g%add_edge(vertices(0),vertices(2),e_ipar,e_rpar)
  edges(3) = g%add_edge(vertices(1),vertices(2),e_ipar,e_rpar)
  edges(4) = g%add_edge(vertices(1),vertices(3),e_ipar,e_rpar)
  edges(5) = g%add_edge(vertices(2),vertices(3),e_ipar,e_rpar)
  edges(6) = g%add_edge(vertices(3),vertices(4),e_ipar,e_rpar)
  edges(7) = g%add_edge(vertices(4),vertices(5),e_ipar,e_rpar)

  ! undirected
  expected_vb = real([0.0,1.5,1.5,6.0,4.0,0.0],dp) / 10.0_dp
  expected_eb = real([2.5,2.5,1.0,4.5,4.5,8.0,5.0],dp) / 15.0_dp
end subroutine test_graph2


subroutine test_graph3(g, is_directed, expected_vb, expected_eb)
  use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
  use iso_fortran_env, only : dp=>real64, output_unit
  implicit none (type, external)
  type(graph_t), intent(inout) :: g
  real(dp), allocatable, intent(out) :: expected_vb(:), expected_eb(:)
  logical, intent(in) :: is_directed


  type(handle_t), allocatable :: edges(:), vertices(:)
  real(dp) :: v_rpar(NRV_PARS), e_rpar(NRE_PARS)
  integer :: v_ipar(NIV_PARS), e_ipar(NIE_PARS)

  print '("Graph 3  directed =",l1)', is_directed
  call g%initialize(is_directed_graph=is_directed)

  allocate(vertices(0:4))
  v_ipar = 1
  v_rpar = 0.0_dp
  v_rpar(4) = 0.5_dp

  v_rpar(1:3) = real([0.0, 1.0, 0.0], dp)
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([1.0, 1.0, 0.0], dp)
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([1.0, 0.0, 0.0], dp)
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([2.0, 1.0, 0.0], dp)
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(1:3) = real([2.0, 0.0, 0.0], dp)
  vertices(4) = g%add_vertex(v_ipar, v_rpar)

  e_ipar = 10
  e_rpar = 0.0_dp
  e_rpar(1) = 1.0_dp
  allocate(edges(5))
  edges(1) = g%add_edge(vertices(0),vertices(1),e_ipar,e_rpar)
  edges(2) = g%add_edge(vertices(0),vertices(2),e_ipar,e_rpar)
  edges(3) = g%add_edge(vertices(1),vertices(3),e_ipar,e_rpar)
  edges(4) = g%add_edge(vertices(2),vertices(3),e_ipar,e_rpar)
  edges(5) = g%add_edge(vertices(3),vertices(4),e_ipar,e_rpar)

  ! directed
  expected_vb = real([0.0,1.0,1.0,3.0,0.0],dp) / 12.0_dp
  expected_eb = real([2.0,2.0,3.0,3.0,4.0],dp) / 20.0_dp
end subroutine test_graph3