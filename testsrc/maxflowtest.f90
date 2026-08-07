  program maxflow_test
    use iso_fortran_env, only : DP=>real64, output_unit, I1=>int8
    use vtuio_mod, only : vtuio_write, vtuio_read, vtuio_data_t
    use graph_mod, only : graph_t, handle_t, MAXFLOW_DINIC, MAXFLOW_EDMOND_KARP
    use graph_user_mod, only : EPOS_WEIGHT, EPOS_FLOW, VPOS_TYPE, VTUIO_MASK
    use graph_testutils_mod, only : testsample_t, parse_lines
    use parse_mod, only : string_t, read_strings
    use utest_mod, only : utest_t
    implicit none (type, external)

    integer, parameter :: algorithm(*) = [MAXFLOW_EDMOND_KARP, MAXFLOW_DINIC]
    integer :: i, j, k, ia
    real(dp) :: maxflow, mincut, flow_cut, flow_cg, flow_cut_cg
    real(dp), allocatable :: edgeflow_cg(:)
    type(vtuio_data_t) :: vtudata
    type(testsample_t) :: ts
    type(string_t), allocatable :: lines(:)
    type(handle_t), allocatable :: s_list(:), t_list(:)
    character(len=2) :: numstr
    character(len=1) :: typestr
    character(len=:), allocatable :: algorithm_str
    type(utest_t) :: utest

    ! Min-cut tests
    print '("GLOBAL MIN-CUT UNIT TESTS")'
    ts%is_directed_graph = .false.
    lines = read_strings('assets/mincut_sample_graphs.txt')
    k = 1
    i = 1
    do while (i <= size(lines))
      write(numstr,'(i2.2)') k
      call parse_lines(lines, i, ts)
      call ts%g%mincut(EPOS_WEIGHT, mincut, s_list=s_list, t_list=t_list)

      ! Test 1 - mincut value match expected
      call utest%assert(mincut, ts%expected_mincut, &
          numstr//' mincut')

      ! Test 2 - S and T lists are correct
      call utest%assert( &
          s_list(:)%get_index_to_map(ts%g), &
          t_list(:)%get_index_to_map(ts%g), &
          ts%expected_s, &
          ts%expected_t, &
          numstr//' ST-lists')

      k = k + 1
#ifdef DEBUG
      print *
#endif
    end do

    ! Max flow tests
    call vtudata%add_item('capacity',start=EPOS_WEIGHT,iclass=3,ncomp=1,nbytes=8)
    call vtudata%add_item('flow',    start=EPOS_FLOW,iclass=3,ncomp=1,nbytes=8)
    lines = read_strings('assets/maxflow_sample_graphs.txt')
    do ia=1,2

      select case(ia)
      case(1)
        print '("EDMOND-KARP UNIT TESTS")'
        algorithm_str = 'Edm-Karp '
      case(2)
        print '("DINIC UNIT TESTS")'
        algorithm_str = 'Dinic '
      end select
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
          write(numstr,'(i2.2)') k
          ! import graph
          call parse_lines(lines, i, ts)
          ! calculate maxflow
          call ts%g%maxflow( &
              ts%g%vertices(ts%sources(1))%handle, &
              ts%g%vertices(ts%sinks(1))%handle, &
              EPOS_WEIGHT, maxflow, &
              position_mincutlabel=VPOS_TYPE, position_flow=EPOS_FLOW, &
              algorithm = algorithm(ia))

          ! Test: maxflow equals expected value
          call utest%assert(maxflow, ts%expected_maxflow(j), &
            algorithm_str//' '//numstr//typestr//' maxflow')

          ! Calculate conductivity (undirected graphs only)
          if (.not. ts%g%is_directed() .and. ia==1) then
            block
              integer, allocatable :: bclabel(:)

              allocate(bclabel(ts%g%nvertices), source=0)
              bclabel(ts%sources(1)) = 2
              bclabel(ts%sinks(1)) = 1
              call ts%g%conductance(EPOS_WEIGHT, bclabel, 0.0_dp, 100.0_dp, &
                  flow_cg, edge_flow=edgeflow_cg)
            end block
          else
            if (allocated(edgeflow_cg)) deallocate(edgeflow_cg)
          end if

          ! Sum flows across min-cut plane
          block
            integer :: ie, iv(2)
            flow_cut = 0.0_dp
            flow_cut_cg = 0.0_dp
            do ie=1,ts%g%nedges
              iv = ts%g%edges(ie)%vertex_indices(ts%g)
              if (is_mincut_edge( &
                  ts%g%vertices(iv(1))%ipar(VPOS_TYPE), &
                  ts%g%vertices(iv(2))%ipar(VPOS_TYPE))) then
                flow_cut = abs(ts%g%edges(ie)%rpar(EPOS_FLOW)) + flow_cut
                if (allocated(edgeflow_cg)) then
                  flow_cut_cg = abs(edgeflow_cg(ie)) + flow_cut_cg
                end if
              end if
            end do

            ! Test: flow across min-cut plane equals maxflow
            call utest%assert(maxflow, flow_cut, &
              algorithm_str//' '//numstr//typestr//' flow across cut')

            ! Test: flow from CG equals the flow across min-cut plane
            if (allocated(edgeflow_cg)) then
              call utest%assert(flow_cg, flow_cut_cg, &
                  'CG solver '//numstr//typestr)
            end if
          end block

          ! write to file for Paraview inspection
          call vtuio_write( &
              'maxflowsample'//numstr//typestr, &
              ts%g, VTUIO_MASK, vtudata=vtudata)

          k=k+1
#ifdef DEBUG
          print *
#endif
        end do ! next graph
      end do   ! next directed/undirected
    end do     ! next algorithm

    call utest%summarize

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
