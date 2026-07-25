  program atom_test
    use iso_fortran_env, only : DP=>real64, output_unit, I1=>int8
    use vtuio_mod, only : vtuio_write, vtuio_read, vtuio_data_t
    use graph_mod, only : graph_t, handle_t, edge_t
    use graph_user_mod, only : VTUIO_MASK, ESIZE_RPAR, ESIZE_IPAR, VSIZE_IPAR, &
        VSIZE_RPAR, VPOS_X, VPOS_VB, VPOS_RADIUS, VPOS_TYPE, EPOS_TYPE
    implicit none (type, external)

    real(DP) :: time
    integer :: i, lab_count

    type(graph_t) :: g, gnew, gbig
    type(handle_t), allocatable :: atom_handles(:), cone_handles(:)
    type(vtuio_data_t) :: vtudata

    real(dp) :: e_rpar(ESIZE_RPAR)
    integer :: e_ipar(ESIZE_IPAR)

    interface
      pure logical function select_edge(this, e)
        use graph_mod, only : edge_t, graph_t
        implicit none (type, external)
        class(graph_t), intent(in) :: this
        type(edge_t), intent(in) :: e
      end function
      pure logical function select_vertex(this, v)
        use graph_mod, only : vertex_t, graph_t
        implicit none (type, external)
        class(graph_t), intent(in) :: this
        type(vertex_t), intent(in) :: v
      end function
    end interface

    ! Initialize graph
    call g%initialize()
   !call g%initialize(is_directed_graph=.true.)

    ! Testing sample of atoms / cones
    allocate(atom_handles(5))

    call add_vertex(g, atom_handles(1), &
      0.75_DP, [0.0_DP, 0.0_DP, 0.0_DP], 1, [1.0_DP, 1.0_DP, 0.0_DP])
    call add_vertex(g, atom_handles(2), &
      0.25_DP, [1.0_DP,0.0_DP,0.0_DP], 1, [1.0_DP, 2.0_DP, 0.0_DP])
    call add_vertex(g, atom_handles(3), &
      0.10_DP, [0.5707_DP,0.6297_DP,0.0_DP], 1, [2.0_DP, 1.0_DP, 0.0_DP])
    call add_vertex(g, atom_handles(4), &
      0.05_DP, [0.0_DP,0.0_DP,1.0_DP], 2, [0.0_DP, 0.0_DP, 0.1_DP])
    call add_vertex(g, atom_handles(5), &
      0.05_DP, [0.0_DP,0.0_DP,2.0_DP], 2, [0.0_DP, 0.0_DP, 0.1_DP])

    allocate(cone_handles(4))
   !allocate(cone_handles(8))
    e_rpar = 0.0_DP
    e_ipar = 0

    e_ipar(EPOS_TYPE) = 10
    cone_handles(1) = g%add_edge(atom_handles(1), atom_handles(2), e_ipar, e_rpar)
   !cone_handles(5) = g%add_edge(atom_handles(2), atom_handles(1), e_ipar, e_rpar)
    e_ipar(EPOS_TYPE) = 10
    cone_handles(2) = g%add_edge(atom_handles(1), atom_handles(3), e_ipar, e_rpar)
   !cone_handles(6) = g%add_edge(atom_handles(3), atom_handles(1), e_ipar, e_rpar)
    e_ipar(EPOS_TYPE) = 20
    cone_handles(3) = g%add_edge(atom_handles(1), atom_handles(4), e_ipar, e_rpar)
   !cone_handles(7) = g%add_edge(atom_handles(4), atom_handles(1), e_ipar, e_rpar)
    e_ipar(EPOS_TYPE) = 20
    cone_handles(4) = g%add_edge(atom_handles(4), atom_handles(5), e_ipar, e_rpar)
   !cone_handles(8) = g%add_edge(atom_handles(5), atom_handles(4), e_ipar, e_rpar)

    ! Write graph to file
    call vtudata%add_item('velocity',start=VPOS_VB,iclass=2,ncomp=3,nbytes=4)
    call vtuio_write('test', g, vtuio_mask, time=123.0_DP, vtudata=vtudata)

    ! Test removing vertex
   !call g%remove_vertex(atom_handles(4))
   !call g%remove_orphaned_edges()
   !call g%print(output_unit)

    ! Read back from file to a new structure
    call vtuio_read('test', gnew, vtuio_mask, time)

    ! Write back a read copy
    call vtuio_write('test_copy', gnew, vtuio_mask, 456.0_DP)

    ! Compare the two
    print *, 'These dumps should be the same'
    call g%print(output_unit)
    call gnew%print(output_unit)

goto 100
    do i=1, min(gnew%nvertices, 100)
      print *, gnew%vertices(i)%rpar
      print *, gnew%vertices(i)%rpar==g%vertices(i)%rpar
      print *, gnew%vertices(i)%ipar, gnew%vertices(i)%ipar==g%vertices(i)%ipar
      print *
    end do
    do i=1, min(gnew%nedges, 100)
      print *, gnew%edges(i)%ipar, gnew%edges(i)%ipar==g%edges(i)%ipar
    end do
100 continue

    ! label connected components
    print *, 'label con com'
    call g%remove_edge(cone_handles(4))
    call g%connected_components( &
        position_label=1, lab_count=lab_count, eselector=select_edge)
    print *, 'Label connected components ', lab_count
    call g%print(output_unit)


    ! remove some objects
    print *, 'Test to remove somethinh...'
    call g%remove_vertex(atom_handles(4))

    call g%print(output_unit)
    call vtuio_write('test_remove', g, vtuio_mask, 123.0_DP)

   !stop 7
    print *
    print *, '*** big file ***'
    call vtuio_read('big', gbig, vtuio_mask)
    call vtuio_write('big_copy', gbig, vtuio_mask)

  contains
    subroutine add_vertex(g, handle, r, x, vtype, vec3)
      type(graph_t), intent(inout) :: g
      type(handle_t), intent(out) :: handle
      real(dp), intent(in) :: r, x(3)
      integer, intent(in) :: vtype
      real(dp), intent(in) :: vec3(3)

      real(dp) :: v_rpar(VSIZE_RPAR)
      integer :: v_ipar(VSIZE_IPAR)

      v_rpar = 0.0_DP
      v_rpar(VPOS_RADIUS) = r
      v_rpar(VPOS_X:VPOS_X+2) = x
      v_rpar(VPOS_VB:VPOS_VB+2) = vec3
      v_ipar = 0
      v_ipar(VPOS_TYPE) = vtype
      handle = g%add_vertex(v_ipar, v_rpar)
    end subroutine add_vertex

  end program atom_test


  pure logical function select_edge(this, e)
    use graph_mod, only : edge_t, graph_t
    implicit none (type, external)
    class(graph_t), intent(in) :: this
    type(edge_t), intent(in) :: e
    select_edge = e%ipar(1)==1 ! just to avouid "unused argument warning"
    select_edge = .true.
   !select_edge = .false.
  end function


  pure logical function select_vertex(this, v)
    use graph_mod, only : vertex_t, graph_t
    implicit none (type, external)
    class(graph_t), intent(in) :: this
    type(vertex_t), intent(in) :: v
    select_vertex = v%ipar(1)==1 ! just to avouid "unused argument warning"
    select_vertex = .true.
  end function
