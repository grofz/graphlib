  program test_meshgen
    use graph_mod, only : graph_handle_t=>handle_t
    use map_mod, only : VPOS_C, VPOS_X
    use mesh_mod, only : mesh_t
    use vtuio_mod, only : vtuio_write, vtuio_data_t, vtuio_read
    use iso_fortran_env, only : dp=>real64
    implicit none

    type(mesh_t) :: m, m2
    type(vtuio_data_t) :: vtudata
    real(dp) :: p0(3), p1(3), p2(3)
    type(graph_handle_t), allocatable :: b(:)
    integer :: i, boffset(0:4)

    call m%initialize(is_3d=.false.)
    print *, 'npoints_per_cell ', m%npoints_per_cell()
    call vtudata%add_item('conc', VPOS_C, 2, 1, 4)

    p0 = 0.0
    p1 = [10.0, 0.0, 0.0]
    p2 = [0.0, 10.0, 0.0]

    call m%append_rectilinear_mesh(p0,p1,p2, 0.5_dp, 0.05_dp, VPOS_X, b, boffset)
    do i=1, m%ncells
      m%vertices(m%index_from_handle(m%cells(i)%dual_vertex))% &
          rpar(VPOS_C) = real(i)
    end do

    ! remove 01 and 23 ghost nodes
    do i=boffset(0)+1, boffset(1)
      call m%remove_vertex(b(i))
    end do
   !do i=boffset(2)+1, boffset(3)
   !  call m%remove_vertex(b(i))
   !end do


    call vtuio_write('meshgen', m, vtudata=vtudata)
    call vtuio_write('meshgenv', m%graph_t, position_id=VPOS_X, vtudata=vtudata)

    call vtuio_read('meshgen', m2, vtudata=vtudata)

  end program test_meshgen
