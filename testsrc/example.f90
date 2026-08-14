program graph_example
  use iso_fortran_env, only : dp=>real64
  use graph_user_mod, only : VSIZE_RPAR, VSIZE_IPAR, ESIZE_RPAR, ESIZE_IPAR
  use graph_mod, only : handle_t, graph_t
  implicit none (type, external)

  ! to store vertex and edge related data
  integer  :: v_ipar(VSIZE_IPAR), e_ipar(ESIZE_IPAR)
  real(dp) :: v_rpar(VSIZE_RPAR), e_rpar(ESIZE_RPAR)

  type(graph_t) :: g
  type(handle_t), allocatable :: vertices(:)
  integer, parameter :: epos_weight=1, epos_maxflow=2, vpos_label=1

  ! define an undirected graph
  call g%initialize(is_directed_graph=.false.)
  block
    integer :: i
    integer, parameter :: nv=6
    type(handle_t) :: edge
    ! add some vertices
    allocate(vertices(nv))
    do i=1, nv 
      vertices(i) = g%add_vertex(v_ipar, v_rpar)
    end do
    ! add some edges
    edge = g%add_edge(vertices(1),vertices(2),e_ipar,e_rpar)
    edge = g%add_edge(vertices(1),vertices(3),e_ipar,e_rpar)
    edge = g%add_edge(vertices(2),vertices(4),e_ipar,e_rpar)
    edge = g%add_edge(vertices(3),vertices(5),e_ipar,e_rpar)
    edge = g%add_edge(vertices(4),vertices(6),e_ipar,e_rpar)
    edge = g%add_edge(vertices(5),vertices(6),e_ipar,e_rpar)
    ! update edges weight
    g%edges(1:g%nedges)%rpar(epos_weight) = &
        real([1.0, 2.0, 0.5, 1.5, 0.1, 4.0], kind=dp)
  end block

  ! Compute conductance of the network
  ! Edge weight is edge conductance
  block
    real(dp), parameter :: x_low=0.0_dp, x_high=100.0_dp
    real(dp), allocatable :: x(:), edge_flow(:)
    integer, allocatable :: bclabel(:)
    real(dp) :: flow, conductance

    allocate(bclabel(g%nvertices), source=0)
    bclabel(1) = 2            ! x_height
    bclabel(g%nvertices) = 1  ! x_low
    call g%conductance(epos_weight, bclabel, x_low, x_high, flow, &
        xfield=x, edge_flow=edge_flow)
    conductance = flow / (x_high - x_low)
    print '("Flow is ",g0,". Conductance is ",g0,".")', flow, conductance
    print '("Nodes potentials: ",/,*(g0,1x))', x
    print '("Flow over edges: ",/,*(g0,1x))', edge_flow
  end block

  ! Compute maximum flow
  ! Edge weight is capacity
  block
    real(dp) :: maxflow
    call g%maxflow(source=vertices(1), sink=vertices(g%nvertices), &
        position_capacity=epos_weight, flow=maxflow, &
        position_flow=epos_maxflow, position_mincutlabel=vpos_label)
    print '("Maximum flow is ",g0,".")', maxflow
    print '("Flow over edges: ",/,*(g0,1x))', &
        g%edges(1:g%nedges)%rpar(epos_maxflow)
    print '("Vertices classification to S and T sets: ",/,*(i3))', &
        g%vertices(1:g%nvertices)%ipar(vpos_label)
  end block

end program