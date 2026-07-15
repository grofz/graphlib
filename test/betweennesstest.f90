! ----------------------------------------------
! A unit test for grpah_betweenness calculation
! ----------------------------------------------
program betweeness
  use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
  use iso_fortran_env, only : dp=>real64, output_unit
  use vtuio_mod, only : vtuio_read, vtuio_write, vtuio_data_t
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
  procedure(test_graph_ai) :: test_graph1, test_graph2, test_graph3, test_graph4
  type(vtuio_data_t) :: vtudata

  integer, parameter :: pos_cost=1, pos_eb=2, pos_vb=5
  integer, parameter :: mask_for_vtuio(*) = [1, 2, 1, 1]

  call vtudata%add_item('edge_betweenness', pos_eb, 3, 1, 4)
  call vtudata%add_item('vertex_betweenness', pos_vb, 2, 1, 4)

  do i=1,4
    ! Sample graphs
    select case(i)
    case(1)
      call test_graph1(g, .false., expected_vb, expected_eb)
    case(2)
      call test_graph2(g, .false., expected_vb, expected_eb)
    case(3)
      call test_graph3(g, .true., expected_vb, expected_eb)
    case(4)
      call test_graph4(g, .false., expected_vb, expected_eb)
    end select

    if (i==4) then
      call g%betweenness(pos_cost, position_eb=pos_eb, position_vb=pos_vb, is_normalized=.false.)
    else
      call g%betweenness(pos_cost, position_eb=pos_eb, position_vb=pos_vb, is_normalized=.true.)
    end if
    print '("Vertex central betweenness: expected/got")'
    print '(*(g0.6,1x))', expected_vb
    print '(*(g0.6,1x))', g%vertices(1:g%nvertices)%rpar(pos_vb)
    print '("Edge central betweenness: expected/got")'
    print '(*(g0.6,1x))', expected_eb
    print '(*(g0.6,1x))', g%edges(1:g%nedges)%rpar(pos_eb)
    print '("Passed? ",l1)', all(g%vertices(1:g%nvertices)%rpar(pos_vb)==expected_vb)
    print '("Passed? ",l1)', all(g%edges(1:g%nedges)%rpar(pos_eb)==expected_eb)
    print *
    if (i==1) call vtuio_write('aa', g, mask_for_vtuio, vtudata=vtudata)
    if (i==4) call vtuio_write('cc', g, mask_for_vtuio, vtudata=vtudata)
  end do
  stop 1

  block
    call g%initialize()
    call vtuio_read('LM50V04', g, mask_for_vtuio)
    g%edges(1:g%nedges)%rpar(pos_cost) = 1.0_dp
    print *, 'Calculating betweenness....'
    call g%betweenness(pos_cost, position_eb=pos_eb, position_vb=pos_vb, is_normalized=.true.)
    print *, '...ok'
    call vtuio_write('bb.vtu', g, mask_for_vtuio, vtudata=vtudata)
  end block

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

  v_rpar(2:4) = real([0.0, 0.0, 0.0], dp)
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([0.0, 1.0, 0.0], dp)
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 1.0, 0.0], dp)
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 0.0, 0.0], dp)
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([2.0, 2.0, 0.0], dp)
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

  v_rpar(2:4) = real([0.0, 1.0, 0.0], dp)
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 1.0, 0.0], dp)
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([0.0, 0.0, 0.0], dp)
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 0.0, 0.0], dp)
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([2.0, 0.0, 0.0], dp)
  vertices(4) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([2.0, 1.0, 0.0], dp)
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

  v_rpar(2:4) = real([0.0, 1.0, 0.0], dp)
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 1.0, 0.0], dp)
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 0.0, 0.0], dp)
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([2.0, 1.0, 0.0], dp)
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([2.0, 0.0, 0.0], dp)
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


subroutine test_graph4(g, is_directed, expected_vb, expected_eb)
  use graph_mod, only : graph_t, handle_t, NIE_PARS, NRE_PARS, NIV_PARS, NRV_PARS
  use iso_fortran_env, only : dp=>real64, output_unit
  implicit none (type, external)
  type(graph_t), intent(inout) :: g
  real(dp), allocatable, intent(out) :: expected_vb(:), expected_eb(:)
  logical, intent(in) :: is_directed


  type(handle_t), allocatable :: edges(:), vertices(:)
  real(dp) :: v_rpar(NRV_PARS), e_rpar(NRE_PARS)
  integer :: v_ipar(NIV_PARS), e_ipar(NIE_PARS)

  print '("Graph 4  directed =",l1)', is_directed
  call g%initialize(is_directed_graph=is_directed)

  allocate(vertices(0:14))
  v_ipar = 1
  v_rpar = 0.0_dp
  v_rpar(1) = 0.5_dp

  v_rpar(2:4) = real([0.0, 0.0, 0.0], dp)   ! 0
  vertices(0) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-1.0, 1.0, 0.0], dp)  ! 1
  vertices(1) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([0.0, -1.0, 0.0], dp)  ! 2
  vertices(2) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-1.0, -1.0, 0.0], dp) ! 3
  vertices(3) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-0.5, -2.0, 0.0], dp) ! 4
  vertices(4) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, 0.5, 0.0], dp)   ! 5
  vertices(5) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-2.0, 1.5, 0.0], dp)  ! 6
  vertices(6) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-2.0, -1.0, 0.0], dp) ! 7
  vertices(7) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-2.3, -0.1, 0.0], dp) ! 8
  vertices(8) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-1.0, 1.9, 0.0], dp)  ! 9
  vertices(9) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, -1.0, 0.0], dp)  ! 10
  vertices(10) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-0.8, -3.0, 0.0], dp) ! 11
  vertices(11) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-2.2, -2.0, 0.0], dp) ! 12
  vertices(12) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([-3.0, -1.0, 0.0], dp) ! 13
  vertices(13) = g%add_vertex(v_ipar, v_rpar)
  v_rpar(2:4) = real([1.0, -2.0, 0.0], dp)  ! 14
  vertices(14) = g%add_vertex(v_ipar, v_rpar)

  e_ipar = 10
  e_rpar = 0.0_dp
  e_rpar(1) = 1.0_dp
  allocate(edges(14))
  edges(1) = g%add_edge(vertices(1),vertices(6),e_ipar,e_rpar)
  edges(2) = g%add_edge(vertices(1),vertices(9),e_ipar,e_rpar)
  edges(3) = g%add_edge(vertices(0),vertices(1),e_ipar,e_rpar)
  edges(4) = g%add_edge(vertices(0),vertices(5),e_ipar,e_rpar)
  edges(5) = g%add_edge(vertices(0),vertices(2),e_ipar,e_rpar)

  edges(6) = g%add_edge(vertices(2),vertices(10),e_ipar,e_rpar)
  edges(7) = g%add_edge(vertices(2),vertices(14),e_ipar,e_rpar)
  edges(8) = g%add_edge(vertices(2),vertices(3),e_ipar,e_rpar)
  edges(9) = g%add_edge(vertices(2),vertices(4),e_ipar,e_rpar)
  edges(10) = g%add_edge(vertices(4),vertices(11),e_ipar,e_rpar)

  edges(11) = g%add_edge(vertices(3),vertices(7),e_ipar,e_rpar)
  edges(12) = g%add_edge(vertices(7),vertices(8),e_ipar,e_rpar)
  edges(13) = g%add_edge(vertices(7),vertices(12),e_ipar,e_rpar)
  edges(14) = g%add_edge(vertices(7),vertices(13),e_ipar,e_rpar)
  ! directed
  expected_vb = real([43,25,70,40,13,0,0,36,0,0,0,0,0,0,0],dp)
  expected_eb = real([0],dp) / 20.0_dp
end subroutine test_graph4