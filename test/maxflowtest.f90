  program maxflow_test
    use iso_fortran_env, only : DP=>real64, output_unit, I1=>int8
    use vtuio_mod, only : vtuio_write, vtuio_read, vtuio_data_t
    use graph_mod, only : graph_t, handle_t
    use graph_user_mod, only : EPOS_WEIGHT
    use graph_testutils_mod, only : testsample_t, parse_lines
    use parse_mod, only : string_t, read_strings
    implicit none (type, external)

    real(DP) :: time
    integer :: i, lab_count, a

    type(graph_t) :: g
    type(handle_t), allocatable :: atom_handles(:), cone_handles(:)
    type(vtuio_data_t) :: vtudata

    integer, parameter :: mask_for_vtuio(*) = [1, 2, 1, 1]
    real(dp) :: maxflow, mincut

    type(testsample_t) :: ts
    type(string_t), allocatable :: lines(:)
    type(handle_t), allocatable :: s_list(:), t_list(:)

    ! Read-file
    lines = read_strings('assets/mincut_sample_graphs.txt')
    i = 1
    do while (i <= size(lines))
      call parse_lines(lines, i, ts)
      call ts%g%mincut(EPOS_WEIGHT, mincut, s_list=s_list, t_list=t_list)
      print '("MINCUT RESULT: ",g0," (expected ",g0,")")', mincut, ts%expected_mincut
      print '("S-LIST ",*(i0,1x))', s_list(:)%get_index_to_map()
      print '("S-LIST EXPECTED ",*(i0,1x))', ts%expected_s
      print '("T-LIST ",*(i0,1x))', t_list(:)%get_index_to_map()
      print '("T-LIST EXPECTED ",*(i0,1x))', ts%expected_t
      print *
    end do
    stop 1


    ! Initialize graph
!   call g%initialize(is_directed_graph=.true.)
    call g%initialize()
!   allocate(atom_handles(4), cone_handles(4))
    allocate(atom_handles(6), cone_handles(9))

    a=1
   !a=2
    select case(a)
    case(1)
      atom_handles(1)=g%add_vertex([1, 0, 0], real([0.75, 1.0,1.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(2)=g%add_vertex([2, 0, 0], real([0.50, 2.0,2.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(3)=g%add_vertex([3, 0, 0], real([0.40, 2.0,0.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(4)=g%add_vertex([4, 0, 0], real([0.90, 3.0,1.0,0.0, 0.0,0.0,0.0],DP))

      cone_handles(1) = g%add_edge(atom_handles(1), atom_handles(2), [10], real([10.0, 0.0],DP))
      cone_handles(2) = g%add_edge(atom_handles(2), atom_handles(3), [20], real([ 5.0, 0.0],DP))
      cone_handles(3) = g%add_edge(atom_handles(3), atom_handles(4), [30], real([11.0, 0.0],DP))
      cone_handles(4) = g%add_edge(atom_handles(1), atom_handles(3), [40], real([ 2.0, 0.0],DP))

      call g%maxflow(atom_handles(1),atom_handles(4),1,maxflow,position_mincutlabel=1)
    case(2)
      atom_handles(1)=g%add_vertex([1, 0, 0], real([0.50, 0.0,1.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(2)=g%add_vertex([1, 0, 0], real([0.75, 1.0,2.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(3)=g%add_vertex([1, 0, 0], real([0.75, 1.0,0.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(4)=g%add_vertex([1, 0, 0], real([0.75, 2.0,2.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(5)=g%add_vertex([1, 0, 0], real([0.75, 2.0,0.0,0.0, 0.0,0.0,0.0],DP))
      atom_handles(6)=g%add_vertex([1, 0, 0], real([0.25, 3.0,1.0,0.0, 0.0,0.0,0.0],DP))

      cone_handles(1) = g%add_edge(atom_handles(1), atom_handles(2), [10], real([10.0, 0.0],DP))
      cone_handles(2) = g%add_edge(atom_handles(1), atom_handles(3), [10], real([10.0, 0.0],DP))
      cone_handles(3) = g%add_edge(atom_handles(2), atom_handles(3), [ 2], real([ 2.0, 0.0],DP))
      cone_handles(4) = g%add_edge(atom_handles(2), atom_handles(4), [ 4], real([ 4.0, 0.0],DP))
      cone_handles(5) = g%add_edge(atom_handles(2), atom_handles(5), [ 8], real([ 8.0, 0.0],DP))
      cone_handles(6) = g%add_edge(atom_handles(3), atom_handles(5), [ 9], real([ 9.0, 0.0],DP))
      cone_handles(7) = g%add_edge(atom_handles(4), atom_handles(5), [ 6], real([ 6.0, 0.0],DP))
      cone_handles(8) = g%add_edge(atom_handles(4), atom_handles(6), [10], real([10.0, 0.0],DP))
      cone_handles(9) = g%add_edge(atom_handles(5), atom_handles(6), [10], real([10.0, 0.0],DP))

      call g%maxflow(atom_handles(1),atom_handles(6),1,maxflow,position_mincutlabel=1,position_flow=2)
      call g%print(output_unit)
    end select
    print *, 'max flow is ',maxflow

    print *, 'STOER-WAGNER'
    call g%mincut(1, mincut, s_list=s_list, t_list=t_list)
    print *, 'min cut is ', mincut



    call vtudata%add_item('capacity',start=1,iclass=3,ncomp=1,nbytes=8)
    call vtudata%add_item('flow',    start=2,iclass=3,ncomp=1,nbytes=8)
    call vtuio_write('aa', g, mask_for_vtuio, time=123.0_DP, vtudata=vtudata)
    call g%print(output_unit)

    ! label connected components
   !print *, 'label con com'
   !call g%remove_edge(cone_handles(4))
   !call g%labconcom(1, lab_count=lab_count, open_edge_f=select_edge)
   !print *, 'Label connected components ', lab_count
   !call g%print(output_unit)

  end program maxflow_test

