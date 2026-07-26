  program maxflow_test
    use iso_fortran_env, only : DP=>real64, output_unit, I1=>int8
    use vtuio_mod, only : vtuio_write, vtuio_read, vtuio_data_t
    use graph_mod, only : graph_t, handle_t
    use graph_user_mod, only : EPOS_WEIGHT, EPOS_FLOW, VPOS_TYPE, VTUIO_MASK
    use graph_testutils_mod, only : testsample_t, parse_lines
    use parse_mod, only : string_t, read_strings
    implicit none (type, external)

    integer :: i, j, k
    real(dp) :: maxflow, mincut, flow_across_cut
    type(vtuio_data_t) :: vtudata
    type(testsample_t) :: ts
    type(string_t), allocatable :: lines(:)
    type(handle_t), allocatable :: s_list(:), t_list(:)
    character(len=2) :: numstr
    character(len=1) :: typestr

    ! Min-cut tests
    print '("GLOBAL MIN-CUT UNIT TESTS")'
    ts%is_directed_graph = .false.
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

    ! Max flow tests
    print '(/,/"EDMOND-KARP UNIT TESTS")'
    call vtudata%add_item('capacity',start=EPOS_WEIGHT,iclass=3,ncomp=1,nbytes=8)
    call vtudata%add_item('flow',    start=EPOS_FLOW,iclass=3,ncomp=1,nbytes=8)
    lines = read_strings('assets/maxflow_sample_graphs.txt')
    do j=1,2 ! directed and undirected
      select case(j)
      case(1)
        ts%is_directed_graph = .true.
        typestr = 'D'
      case(2)
        ts%is_directed_graph = .false.
        typestr = 'U'
      end select
      i = 1
      k = 1
      do while (i <= size(lines)) ! for all graphs in the file
        ! import graph
        call parse_lines(lines, i, ts)
        ! calculate maxflow (Edmond-Carp)
        call ts%g%maxflow( &
            ts%g%vertices(ts%sources(1))%handle, &
            ts%g%vertices(ts%sinks(1))%handle, &
            EPOS_WEIGHT, maxflow, &
            position_mincutlabel=VPOS_TYPE, position_flow=EPOS_FLOW)
        print '("MAXFLOW RESULT: ",g0," (expected ",g0,")")', maxflow, ts%expected_maxflow(j)
        ! sum flow across min-cut plane
        block
          integer :: ie, iv(2)
          flow_across_cut = 0.0_dp
          do ie=1,ts%g%nedges
            iv = ts%g%edges(ie)%vertex_indices(ts%g)
            if (is_mincut_edge( &
                ts%g%vertices(iv(1))%ipar(VPOS_TYPE), &
                ts%g%vertices(iv(2))%ipar(VPOS_TYPE))) then
              flow_across_cut = abs(ts%g%edges(ie)%rpar(EPOS_FLOW)) + flow_across_cut
            end if
          end do
          print '("flow across cut ",g0)',flow_across_cut
          print '("PASSED? ",l2)', flow_across_cut==maxflow .and. maxflow==ts%expected_maxflow(j) 
        end block
        ! Write to file for Paraview inspection
        write(numstr,'(i2.2)') k
        call vtuio_write( &
            'maxflowsample'//numstr//typestr, ts%g, VTUIO_MASK, vtudata=vtudata)
        k=k+1
        print *
      end do ! next graph
    end do   ! next directed/undirected
  
  contains
    logical function is_mincut_edge(a, b)
      integer, intent(in) :: a, b

      is_mincut_edge = .false.
      if (a/=1 .and. a/=2) return
      if (b/=1 .and. b/=2) return
      if (a==b) return
      is_mincut_edge = .true.
    end function

  end program maxflow_test

