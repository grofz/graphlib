program scc
  use graph_testutils_mod, only : testsample_t, parse_lines
  use parse_mod, only : string_t, read_strings
  use graph_user_mod, only : VSIZE_IPAR, VSIZE_RPAR, ESIZE_IPAR, ESIZE_RPAR, &
      VTUIO_MASK, VPOS_TYPE
  use graph_mod, only : graph_t, handle_t
  use iso_fortran_env, only : dp=>real64, output_unit
  use vtuio_mod, only : vtuio_read, vtuio_write, vtuio_data_t
  implicit none (type, external)

  type(testsample_t) :: ts
  type(string_t), allocatable :: lines(:)
  integer :: i, k
  character(len=2) :: numstr
  logical, allocatable :: vmask(:), emask(:)
  integer, allocatable :: labels(:), levels(:)
  type(handle_t) :: handle
  integer :: vipar(VSIZE_IPAR)
  real(dp) :: vrpar(VSIZE_RPAR)

  ! Empty graph test
  print '("# Empty graph")'
  call ts%g%initialize(is_directed_graph=.true.)
  call ts%g%strongly_connected_components(labels=labels)
  print '("PASSED ",l2)', size(labels)==0
  call ts%g%topological_levels(levels, components=labels)
  print '("LEVELS PASSED ",l2,/)', ts%g%verify_topological_levels(levels, components=labels)

  ! Single isolated vertex
  print '("# Single isolated vertex")'
  call ts%g%initialize(is_directed_graph=.true.)
  handle = ts%g%add_vertex(vipar, vrpar)
  call ts%g%strongly_connected_components(labels=labels)
  print '("PASSED ",l2)', size(labels)==1 .and. labels(1)==1
  call ts%g%topological_levels(levels, components=labels)
  print '("LEVELS PASSED ",l2,/)', ts%g%verify_topological_levels(levels, components=labels)

  lines = read_strings('assets/scc_sample_graphs.txt')
  i = 1
  k = 1

  do while (i<= size(lines))
    ts%is_directed_graph = .true.
    call parse_lines(lines, i, ts)
    call ts%g%strongly_connected_components(labels=labels )
  print *, size(labels), size(ts%expected_scc)
    print '("PASSED ",l2)', all(labels==ts%expected_scc)
    if (.not. all(labels==ts%expected_scc)) then
      print '("Lables: ",*(i0,1x))', labels
      print '("Expect: ",*(i0,1x))', ts%expected_scc
    end if

    ! Test levels
    call ts%g%topological_levels(levels, components=labels)
    print '("LEVELS PASSED ",l2,/)', ts%g%verify_topological_levels(levels, components=labels)

    write(numstr,'(i2.2)') k
!   call vtuio_write('sccsample'//numstr, ts%g, VTUIO_MASK, vtudata=vtudata)
    k = k+1

    print *
  end do





end program


