program test_handle
  use graph_mod, only : graph_t, handle_t, MAP_NULL
  use mesh_mod, only : mesh_t
  use utest_mod, only : utest_t
  use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR
  use iso_fortran_env, only : dp=>real64, I1B=>int8
  implicit none

  type(utest_t) :: utest
  type(graph_t) :: g1, g2
  type(mesh_t) :: m1, m2
  type(handle_t) :: h1, h2, h1new, j1, j2, hnull
  integer :: ipar(VSIZE_IPAR)
  real(dp) :: rpar(VSIZE_RPAR), xpos(3)

  ! graph_t test
  call g1%initialize()
  call g2%initialize()

  h1 = g1%add_vertex(ipar, rpar)
  h2 = g2%add_vertex(ipar, rpar)

  call utest%assert(g1%index_from_handle(h1)==MAP_NULL, .false., 'graph1 - handle 1 accepted')
  call utest%assert(g2%index_from_handle(h2)==MAP_NULL, .false., 'graph2 - handle 2 accepted')
  call utest%assert(g1%index_from_handle(h2)==MAP_NULL, .true., 'graph1 - handle 2 rejected')
  call utest%assert(g2%index_from_handle(h1)==MAP_NULL, .true., 'graph2 - handle 1 rejected')
  call g1%initialize()
  h1new = g1%add_vertex(ipar, rpar)
  call utest%assert(g1%index_from_handle(h1)==MAP_NULL, .true., 'graph1 - old handle rejected')
  call utest%assert(g1%index_from_handle(h1new)/=MAP_NULL, .true., 'graph1 - new handle accepted')

  ! mesh_t test
  call m1%initialize()
  call m2%initialize()
  xpos = 0.0_dp

  j1 = m1%add_point(xpos)
  j2 = m2%add_point(xpos)

  call utest%assert(m1%index_from_handle(j1)==MAP_NULL, .false., 'mesh1 - handle 1 accepted')
  call utest%assert(m2%index_from_handle(j2)==MAP_NULL, .false., 'mesh2 - handle 2 accepted')
  call utest%assert(m1%index_from_handle(j2)==MAP_NULL, .true., 'mesh1 - handle 2 rejected')
  call utest%assert(m2%index_from_handle(j1)==MAP_NULL, .true., 'mesh2 - handle 1 rejected')
  call utest%assert(m1%index_from_handle(h1)==MAP_NULL, .true., 'mesh1 - graph handle rejected')
  call utest%assert(m1%index_from_handle(h1new)==MAP_NULL, .true., 'mesh1 - new graph handle rejected')
  call utest%assert(m2%index_from_handle(h2)==MAP_NULL, .true., 'mesh2 - graph handle rejected')

  ! null handle
  hnull = handle_t(id=1, type=2_I1B, version=1)
  call utest%assert(g1%index_from_handle(hnull)==MAP_NULL, .true., 'graph1 - null handle rejected')
  call utest%assert(m1%index_from_handle(hnull)==MAP_NULL, .true., 'mesh1 - null handle rejected')

  call utest%summarize()
  if (.not. utest%all_passed()) stop 2
end program