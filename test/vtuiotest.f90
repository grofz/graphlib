  program atom_test
    use iso_fortran_env, only : DP=>real64
    use vtuio_mod, only : vtuio_write, vtuio_read
!   use import_dem2011_mod, only : import_dem2011
    use graph_mod, only : graph_t, handle_t, graph_add_vertex, graph_add_edge
    implicit none (type, external)

    real(DP) :: velo(3), time
    real(DP), allocatable :: x(:,:), x2(:,:)
    integer :: i

    type(graph_t) :: g, gnew, gbig
    type(handle_t), allocatable :: atom_handles(:), cone_handles(:)

    integer, parameter :: mask_for_vtuio(*) = [1, 2, 1, 1]

    ! Initialize graph
    call g%initialize()

    ! Testing sample of atoms / cones
    allocate(atom_handles(5), x(3,5))
    velo = [1.0, 1.0, 0.0]
    x(:,1) = real([0.0,0.0,0.0],DP)
    atom_handles(1) = graph_add_vertex(g, [1], [0.75_DP, x(:,1)])

    velo = [1.0, 2.0, 0.0]
    x(:,2) = real([1.0,0.0,0.0],DP)
    atom_handles(2) = graph_add_vertex(g, [1], [0.25_DP, x(:,2)])

    velo = [2.0, 1.0, 0.0]
    x(:,3) = real([0.5707,0.6297,0.0],DP)
    atom_handles(3) = graph_add_vertex(g, [1], [0.10_DP, x(:,3)])

    velo = [0.0, 0.0, 0.1]
    x(:,4) = real([0.0,0.0,1.0],DP)
    atom_handles(4) = graph_add_vertex(g, [2], [0.05_DP, x(:,4)])

    x(:,5) = real([0.0,0.0,2.0],DP)
    atom_handles(5) = graph_add_vertex(g, [2], [0.05_DP, x(:,5)])

    allocate(cone_handles(4))
    cone_handles(1) = graph_add_edge(g, atom_handles(1), atom_handles(2), [10], [real(dp)::])
    cone_handles(2) = graph_add_edge(g, atom_handles(1), atom_handles(3), [10], [real(dp)::])
    cone_handles(3) = graph_add_edge(g, atom_handles(1), atom_handles(4), [20], [real(dp)::])
    cone_handles(4) = graph_add_edge(g, atom_handles(4), atom_handles(5), [20], [real(dp)::])

    call vtuio_write('test', g, mask_for_vtuio, 123.0_DP)
    print *, 'atom_test write finished'

    call vtuio_read('test', gnew, mask_for_vtuio, time)
    print *, 'atom_test read finished'

    do i=1, min(gnew%nvertices, 100)
      print *, gnew%vertices(i)%rpar(2:4)
      print *, gnew%vertices(i)%rpar(2:4) == g%vertices(i)%rpar(2:4)
      print *, gnew%vertices(i)%rpar(1), g%vertices(i)%rpar(1), gnew%vertices(i)%rpar(1) == g%vertices(i)%rpar(1)
      print *, gnew%vertices(i)%ipar(1), gnew%vertices(i)%ipar(1) == g%vertices(i)%ipar(1)
      print *
    end do
    do i=1, min(gnew%nedges, 100)
      print *, gnew%edges(i)%ipar(1), gnew%edges(i)%ipar(1) == g%edges(i)%ipar(1)
    end do
    call vtuio_write('test_copy', gnew, mask_for_vtuio, 456.0_DP)

    print *
    print *, '*** big file ***'
    call vtuio_read('big', gbig, mask_for_vtuio)
    call vtuio_write('big_copy', gbig, mask_for_vtuio)

!   g = import_dem2011('tmp')
!   call g%writevtu('tmpvtu', rem_dbl_edges=.true.)

  end program atom_test
