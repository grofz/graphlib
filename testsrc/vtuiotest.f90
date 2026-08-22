  program vtuio_test
    use iso_fortran_env, only : SP=>real32, DP=>real64, output_unit, I1=>int8
    use vtuio_mod, only : vtuio_write, vtuio_read, vtuio_data_t, META_IS_CELL, &
      META_IS_POINT, META_IS_INT, META_IS_REAL
    use graph_mod, only : graph_t, graph_handle_t=>handle_t, edge_t
    use mesh_mod, only : mesh_t, mesh_handle_t
    use graph_user_mod, only : ESIZE_RPAR, ESIZE_IPAR, VSIZE_IPAR, VSIZE_RPAR
    use map_mod
    use utest_mod, only : utest_t
    implicit none (type, external)

    type(graph_t) :: gold, gnew, gbig
    type(mesh_t) :: mold, mnew
    type(graph_handle_t), allocatable :: v_handles(:), e_handles(:)
    type(mesh_handle_t), allocatable :: p_handles(:), c_handles(:)
    type(vtuio_data_t) :: vtudata
    type(utest_t) :: utest
    real(SP) :: time
    real(SP), parameter :: time0 = 123.45
    real(dp), allocatable :: vpos(:,:), vrad(:), vvelo(:,:), eflow(:)
    integer, allocatable :: vtype(:), etype(:)

    ! Initialize graph
    call gold%initialize()

    ! Construct graph
    block
      real(dp) :: e_rpar(ESIZE_RPAR), v_rpar(VSIZE_RPAR)
      integer :: e_ipar(ESIZE_IPAR), v_ipar(VSIZE_IPAR)
      integer :: i
      integer, allocatable :: connections(:,:)

      allocate(v_handles(5), e_handles(4))

      ! graph attributes
      vtype = [1, 1, 1, 2, 2] 
      vrad = [0.75_DP, 0.25_DP, 0.10_DP, 0.05_DP, 0.05_DP]
      vpos = reshape( [                 &
          1.0_DP,    0.0_DP,    0.0_DP, &
          0.5707_DP, 0.6297_DP, 0.0_DP, &
          0.0_DP,    0.0_DP,    1.0_DP, &
          0.0_DP,    0.0_DP,    2.0_DP, &
          0.0_DP,    0.0_DP,    0.0_DP  &
          ], [3, size(v_handles)])
      vvelo = reshape( [ &
          1.0_DP, 1.0_DP, 0.0_DP, &
          1.0_DP, 2.0_DP, 0.0_DP, &
          2.0_DP, 1.0_DP, 0.0_DP, &
          0.0_DP, 0.0_DP, 0.1_DP, &
          0.0_DP, 0.0_DP, 0.1_DP  &
          ], [3, size(v_handles)])
      etype = [10, 20, 30, 40]
      eflow = [11.0, 22.0, 33.0, 44.0]
      connections = reshape( [ &
          1, 2, &
          1, 3, &
          1, 4, &
          4, 5  &
          ], [2, size(e_handles)])

      ! graph vertices
      do i=1, size(v_handles)
        v_ipar(VPOS_TYPE)   = vtype(i)
        v_rpar(VPOS_RADIUS) = vrad(i)
        v_rpar(VPOS_X:VPOS_X+2) = vpos(:,i)
        v_rpar(VPOS_VB:VPOS_VB+2) = vvelo(:,i)
        v_handles(i) = gold%add_vertex(v_ipar,v_rpar)
      end do

      ! graph edges
      do i=1, size(e_handles)
        e_ipar(EPOS_TYPE) = etype(i)
        e_rpar(EPOS_FLOW) = eflow(i)
        e_handles(i) = gold%add_edge(v_handles(connections(1,i)), &
            v_handles(connections(2,i)), e_ipar, e_rpar)
      end do
    end block

    ! Write graph
    call vtudata%add_item('velocity', VPOS_VB, META_IS_POINT+META_IS_REAL, 3, 4)
    call vtudata%add_item('radius', VPOS_RADIUS, META_IS_POINT+META_IS_REAL, 1, 8)
    call vtudata%add_item('ctype', EPOS_TYPE, META_IS_CELL+META_IS_INT, 1, 1)
    call vtudata%add_item('vtype', VPOS_TYPE, META_IS_POINT+META_IS_INT, 1, 1)
    call vtudata%add_item('flow', EPOS_FLOW, META_IS_CELL+META_IS_REAL, 1, 8)
    call vtuio_write('test_graph', gold, position_id=VPOS_X, time=real(time0,dp), vtudata=vtudata)

    ! Read graph
    block
      real(dp) :: timeDP
      call vtuio_read('test_graph', gnew, position_id=VPOS_X, time=timeDP, vtudata=vtudata)
      time = real(timeDP)
    end block

    ! Compare gold and gnew
    block
      call utest%assert(time0, time, 'writing / reading time to VTU file works')
      call utest%assert(gold%nvertices, gnew%nvertices, 'number of vertices match')
      call utest%assert(gold%nedges, gnew%nedges, 'number of edges match')
      call utest%assert(.true., &
          all(real(gold%vertices(1:gold%nvertices)%rpar(VPOS_X)) == &
              real(gnew%vertices(1:gnew%nvertices)%rpar(VPOS_X))), 'position x match')
      call utest%assert(.true., &
          all(real(gold%vertices(1:gold%nvertices)%rpar(VPOS_X+1)) == &
              real(gnew%vertices(1:gnew%nvertices)%rpar(VPOS_X+1))), 'position y match')
      call utest%assert(.true., &
          all(real(gold%vertices(1:gold%nvertices)%rpar(VPOS_X+2)) == &
              real(gnew%vertices(1:gnew%nvertices)%rpar(VPOS_X+2))), 'position z match')
      call utest%assert(.true., &
          all((gold%vertices(1:gold%nvertices)%rpar(VPOS_RADIUS)) == &
              (gnew%vertices(1:gnew%nvertices)%rpar(VPOS_RADIUS))), 'radius match')
      call utest%assert(.true., &
          all(real(gold%vertices(1:gold%nvertices)%rpar(VPOS_VB)) == &
              real(gnew%vertices(1:gnew%nvertices)%rpar(VPOS_VB))), 'velocity x match')
      call utest%assert(.true., &
          all(real(gold%vertices(1:gold%nvertices)%rpar(VPOS_VB+1)) == &
              real(gnew%vertices(1:gnew%nvertices)%rpar(VPOS_VB+1))), 'velocity y match')
      call utest%assert(.true., &
          all(real(gold%vertices(1:gold%nvertices)%rpar(VPOS_VB+2)) == &
              real(gnew%vertices(1:gnew%nvertices)%rpar(VPOS_VB+2))), 'velocity z match')
      call utest%assert(.true., &
          all(gold%vertices(1:gold%nvertices)%ipar(VPOS_TYPE) == &
              gnew%vertices(1:gnew%nvertices)%ipar(VPOS_TYPE)), 'vertex type match')
      call utest%assert(.true., &
          all(gold%edges(1:gold%nedges)%ipar(EPOS_TYPE) == &
              gnew%edges(1:gnew%nedges)%ipar(EPOS_TYPE)), 'edges type match')
      call utest%assert(.true., &
          all(gold%edges(1:gold%nedges)%rpar(EPOS_FLOW) == &
              gnew%edges(1:gnew%nedges)%rpar(EPOS_FLOW)), 'flow match')
    end block
#ifdef DEBUG
    call gold%print(output_unit)
    call gnew%print(output_unit)
#endif

    ! test reading a big file
    print *
    print *, '*** big file ***'
    call vtuio_read('big', gbig, position_id=VPOS_X, vtudata=vtudata)
    call vtuio_write('big_copy', gbig, position_id=VPOS_X)

    ! test summary
    call utest%summarize()
    if (.not. utest%all_passed()) stop 1

  end program vtuio_test