  program test_meshgen
    use graph_mod, only : graph_handle_t=>handle_t
    use map_mod, only : VPOS_C, VPOS_X, VPOS_RTMP, EPOS_RTMP, VPOS_BC
    use mesh_mod, only : mesh_t, integrate_pde
    use vtuio_mod, only : vtuio_write, vtuio_data_t, vtuio_read
    use iso_fortran_env, only : dp=>real64
    implicit none

    type(mesh_t) :: m, m2
    type(vtuio_data_t) :: vtudata
    real(dp) :: p0(3), p1(3), p2(3)
    type(graph_handle_t), allocatable :: b(:)
    integer :: i, j, boffset(0:4)
    real(dp), allocatable :: u_init(:), u_out(:,:)
    integer, allocatable :: bc_label(:)
    real(dp) :: t_start, t_end, dt_out, dt_comp

    call m%initialize(is_3d=.false.)
    print *, 'npoints_per_cell ', m%npoints_per_cell()
    call vtudata%add_item('conc', VPOS_C, 2, 1, 4)


    block
      real(dp) :: rel_shift = 0.0_dp
      real(dp) :: cell_size = 0.5_dp

      p0 = 0.0
      p1 = [10.0, 0.0, 0.0]
      p2 = [0.0, 10.0, 0.0]

      call m%append_rectilinear_mesh(p0, p1, p2, cell_size, rel_shift, VPOS_X, &
          b, boffset)
      do i=1, m%ncells
        m%vertices(m%index_from_handle(m%cells(i)%dual_vertex))% &
            rpar(VPOS_C) = real(i)
      end do
    end block

    ! remove 01 and 23 ghost nodes
    do i=boffset(0)+1, boffset(1)
      call m%remove_vertex(b(i))
    end do

    ! set bc_label
    allocate(bc_label(m%nvertices), source=0)
    do i=boffset(2)+1, boffset(3)
      bc_label(m%index_from_handle(b(i))) = 2
    end do

    do i=boffset(1)+1, boffset(2)
      bc_label(m%index_from_handle(b(i))) = 1
    end do
    do i=boffset(3)+1, boffset(4)
      bc_label(m%index_from_handle(b(i))) = 2
    end do

    allocate(u_init(m%nvertices))
    where(bc_label==2)
      u_init = 100.0_dp
    else where (bc_label==1)
      u_init = 0.0_dp
    else where
      u_init = 0.0_dp
    end where

    ! set conductance
    m%edges(1:m%nedges)%rpar(EPOS_RTMP) = 1.0_dp
    ! set capacitance
    m%vertices(1:m%nvertices)%rpar(VPOS_RTMP) = 1.0_dp

    ! integrate
    t_start = 0.0
    t_end = 10.0
    dt_comp = 0.05
    dt_out = 0.5
    call integrate_pde(m, t_start, t_end, dt_comp, dt_out, EPOS_RTMP, VPOS_RTMP, &
      bc_label, u_init, u_out)

    call vtuio_write('meshgen', m, vtudata=vtudata)
    call vtuio_write('meshgenv', m%graph_t, position_id=VPOS_X, vtudata=vtudata)

    ! write uout
    block
      character(len=3) :: text
      do i=1, size(u_out,2)
        write(text,'(i3.3)') i
        m%vertices(1:m%nvertices)%rpar(VPOS_C) = u_out(:,i)
        call vtuio_write('time'//text, m, vtudata=vtudata)
print *, minval(u_out(:,i),1), minloc(u_out(:,i),1), maxval(u_out(:,i),1), maxloc(u_out(:,i),1)
      end do
    end block
!   call vtuio_write('meshgenv', m%graph_t, position_id=VPOS_X, vtudata=vtudata)

    call vtuio_read('meshgen', m2, vtudata=vtudata)

  end program test_meshgen
