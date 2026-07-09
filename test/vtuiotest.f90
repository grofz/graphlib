  program atom_test
    use iso_fortran_env, only : DP=>real64, output_unit
    use vtuio_mod, only : vtuio_write, vtuio_read
!   use import_dem2011_mod, only : import_dem2011
    use graph_mod, only : graph_t, handle_t
    implicit none (type, external)

    real(DP) :: velo(3), time
    real(DP), allocatable :: x(:,:), x2(:,:)
    integer :: i, lab_count

    type(graph_t) :: g, gnew, gbig
    type(handle_t), allocatable :: atom_handles(:), cone_handles(:)

    integer, parameter :: mask_for_vtuio(*) = [1, 2, 1, 1]

    interface
      subroutine graph_export(graph)
        use graph_mod, only : graph_t, handle_t
        implicit none (type, external)
        type(graph_t), intent(in) :: graph
      end subroutine

      logical function select_edge(e)
        use graph_mod, only : edge_t
        implicit none (type, external)
        type(edge_t), intent(in) :: e
      end function
      logical function select_vertex(v)
        use graph_mod, only : vertex_t
        implicit none (type, external)
        type(vertex_t), intent(in) :: v
      end function
    end interface

    ! Initialize graph
    call g%initialize()

    ! Testing sample of atoms / cones
    allocate(atom_handles(5), x(3,5))
    velo = [1.0, 1.0, 0.0]
    x(:,1) = real([0.0,0.0,0.0],DP)
    atom_handles(1) = g%add_vertex([1], [0.75_DP, x(:,1)])

    velo = [1.0, 2.0, 0.0]
    x(:,2) = real([1.0,0.0,0.0],DP)
    atom_handles(2) = g%add_vertex([1], [0.25_DP, x(:,2)])

    velo = [2.0, 1.0, 0.0]
    x(:,3) = real([0.5707,0.6297,0.0],DP)
    atom_handles(3) = g%add_vertex([1], [0.10_DP, x(:,3)])

    velo = [0.0, 0.0, 0.1]
    x(:,4) = real([0.0,0.0,1.0],DP)
    atom_handles(4) = g%add_vertex([2], [0.05_DP, x(:,4)])

    x(:,5) = real([0.0,0.0,2.0],DP)
    atom_handles(5) = g%add_vertex([2], [0.05_DP, x(:,5)])

    allocate(cone_handles(4))
    cone_handles(1) = g%add_edge(atom_handles(1), atom_handles(2), [10], [real(dp)::])
    cone_handles(2) = g%add_edge(atom_handles(1), atom_handles(3), [10], [real(dp)::])
    cone_handles(3) = g%add_edge(atom_handles(1), atom_handles(4), [20], [real(dp)::])
    cone_handles(4) = g%add_edge(atom_handles(4), atom_handles(5), [20], [real(dp)::])

    call vtuio_write('test', g, mask_for_vtuio, 123.0_DP)


    ! Read back from file to a new structure
    call vtuio_read('test', gnew, mask_for_vtuio, time)

    ! Write back a read copy
    call vtuio_write('test_copy', gnew, mask_for_vtuio, 456.0_DP)

    ! Compare the two
    print *, 'These dumps should be the same'
    call g%print(output_unit)
    call gnew%print(output_unit)

goto 100
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
100 continue

    ! label connected components
    print *, 'label con com'
    call g%remove_edge(cone_handles(4))
    call g%labconcom(1, lab_count=lab_count, open_edge_f=select_edge)
    print *, 'Label connected components ', lab_count
    call g%print(output_unit)
    stop 8


    ! remove some objects
    print *, 'Test to remove somethinh...'
    call g%remove_vertex(atom_handles(4))
   !call g%remove_edge(cone_handles(4))

    call g%print(output_unit)
    call vtuio_write('test_remove', g, mask_for_vtuio, 123.0_DP)

!   call graph_export(g)
!   call graph_export(gnew)

   !stop 7
    print *
    print *, '*** big file ***'
    call vtuio_read('big', gbig, mask_for_vtuio)
    call vtuio_write('big_copy', gbig, mask_for_vtuio)

!   g = import_dem2011('tmp')
!   call g%writevtu('tmpvtu', rem_dbl_edges=.true.)
!   call graph_export(gbig)
!   call gbig%print(output_unit)

  end program atom_test


  subroutine graph_export(graph)
    use iso_fortran_env, only : DP=>real64
    use graph_mod, only : graph_t, handle_t
    implicit none (type, external)
    type(graph_t), intent(in) :: graph
    integer :: i
    print *,' VERTICES'
    do i=1, graph%nvertices
      print *, i, graph%vertices(i)%ipar, graph%vertices(i)%rpar
      print *
    end do
    print *,' EDGES'
    do i=1, graph%nedges
      print *, graph%edges(i)%vertex_indices(graph), graph%edges(i)%ipar, graph%edges(i)%rpar
      print *
    end do
  end subroutine graph_export

      logical function select_edge(e)
        use graph_mod, only : edge_t
        implicit none (type, external)
        type(edge_t), intent(in) :: e
        select_edge = .true.
       !select_edge = .false.
      end function
      logical function select_vertex(v)
        use graph_mod, only : vertex_t
        implicit none (type, external)
        type(vertex_t), intent(in) :: v
        select_vertex = .true.
      end function
