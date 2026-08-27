  program test_meshgen
    use map_mod, only : VPOS_C
    use mesh_mod, only : mesh_t
    use vtuio_mod, only : vtuio_write, vtuio_data_t, vtuio_read
    use iso_fortran_env, only : dp=>real64
    implicit none

    type(mesh_t) :: m, m2
    type(vtuio_data_t) :: vtudata
    real(dp) :: p0(3), p1(3), p2(3)
    integer :: i

    call m%initialize(is_3d=.false.)
    print *, 'npoints_per_cell ', m%npoints_per_cell()
    call vtudata%add_item('conc', VPOS_C, 2, 1, 4)

    p0 = 0.0
    p1 = [100.0, 0.0, 0.0]
    p2 = [10.0, 100.0, 0.0]

    call m%append_rectilinear_mesh(p0,p1,p2, 0.5_dp, 0.05_dp)
    do i=1, m%ncells
      m%vertices(m%index_from_handle(m%cells(i)%dual_vertex))% &
          rpar(VPOS_C) = real(i)
    end do
    call vtuio_write('meshgen', m, vtudata=vtudata)

    call vtuio_read('meshgen', m2, vtudata=vtudata)

  end program test_meshgen
