  program test_meshgen
    use graph_mod, only : graph_handle_t=>handle_t
    use map_mod, only : VPOS_C, VPOS_X, VPOS_RTMP, VPOS_VB, EPOS_RTMP, VPOS_BC
    use mesh_mod, only : mesh_t, integrate_pde
    use vtuio_mod, only : vtuio_write, vtuio_data_t, vtuio_read
    use analytical_mod, only : analytical_1d_temperature
    use iso_fortran_env, only : dp=>real64
    implicit none (type, external)

    type(mesh_t) :: m, m2
    type(vtuio_data_t) :: vtudata
    integer :: i, j, fid
    real(dp), allocatable :: u_init(:), u_out(:,:)
    integer, allocatable :: bc_label(:)
    real(dp) :: t_start, t_end, dt_out, dt_comp

    real(dp) :: conductivity, capacity, width, temp_bc(2)

    call m%initialize(is_3d=.false.)
    print *, 'npoints_per_cell ', m%npoints_per_cell()

    ! What will be written to .vtu file for Paraview?
    call vtudata%add_item('conc', VPOS_C, 2, 1, 4)

    ! Parameters
    conductivity = 2.0_dp
    capacity = 1.0_dp
    width = 10.0_dp
    temp_bc(1) = 0.0_dp
    temp_bc(2) = 100.0_dp

    ! Generate mesh
    block
      real(dp) :: rel_shift = 0.0_dp
      real(dp) :: cell_size = 0.1_dp
      real(dp) :: p0(3), p1(3), p2(3)
      integer :: boffset(0:4)
      type(graph_handle_t), allocatable :: b(:)

      p0 = 0.0_dp
      p1 = [width, 0.0_dp, 0.0_dp]
      p2 = [0.0_dp, width, 0.0_dp]

      call m%append_rectilinear_mesh(p0, p1, p2, cell_size, rel_shift, VPOS_X, &
          b, boffset)
      do i=1, m%ncells
        m%vertices(m%index_from_handle(m%cells(i)%dual_vertex))% &
            rpar(VPOS_C) = real(i)
      end do

      ! remove 02 and 13 ghost nodes (west / east)
      do i=boffset(1)+1, boffset(2)
        call m%remove_vertex(b(i))
      end do
      do i=boffset(3)+1, boffset(4)
        call m%remove_vertex(b(i))
      end do

      ! set boundary and initial conditions
      allocate(bc_label(m%nvertices), source=0)
      do i=boffset(2)+1, boffset(3) ! north border
        bc_label(m%index_from_handle(b(i))) = 2 ! temp_2
      end do
      do i=boffset(0)+1, boffset(1) ! south border
        bc_label(m%index_from_handle(b(i))) = 1 ! temp_1
      end do
    end block

    allocate(u_init(m%nvertices))
    where(bc_label==2)
      u_init = temp_bc(2)
    else where (bc_label==1)
      u_init = temp_bc(1)
    else where
      u_init = temp_bc(1)
    end where

    ! set conductivity
    m%vertices(1:m%nvertices)%rpar(VPOS_VB) = conductivity
    ! set capacitt
    m%vertices(1:m%nvertices)%rpar(VPOS_RTMP) = capacity

    ! integrate
    t_start = 0.0
    t_end = 5.0
   !t_end = 10.0
    dt_comp = 0.025
    dt_out = 1.0
    call integrate_pde(m, t_start, t_end, dt_comp, dt_out, &
      VPOS_VB, VPOS_RTMP, EPOS_RTMP, &
      bc_label, u_init, u_out)

    ! compare numerical and analytical solution
    block
      integer :: itime, icount
      real(dp) :: time, positions(3,4), y, u_anal, u_norm, u_err, u_maxerr

      open(newunit=fid, file='convergence.dat', status='replace')
      u_norm = (temp_bc(1) + temp_bc(2)) / 2.0_dp
      do itime=1, size(u_out,2)
        time = t_start + min(t_end, real(itime-1)*dt_out)
        u_err = 0.0
        u_maxerr = 0.0
        icount = 0
        do i=1, m%nvertices
          if (m%v2c_index(i)<1) cycle
          positions = m%cells(m%v2c_index(i))%point_positions(m)
          y = sum(positions(2,1:3)) / 3.0_dp
         !y = m%vertices(i)%rpar(VPOS_X+1)
          u_anal = analytical_1d_temperature(width-y, time, conductivity/capacity, width, temp_bc)
!print '("vid ",i0,", y",g0,", T ",g0,", T_an ",g0)', &
!  i,y,u_out(i,itime), u_anal
          u_err = u_err + abs(u_anal - u_out(i,itime))
          u_maxerr = max(u_maxerr, abs(u_anal - u_out(i,itime)))
          icount = icount + 1
          write(fid,*) time, y, u_anal, u_out(i,itime)
        end do
        write(fid,*)
        print '("time = ",g0,"   avg.error = ",e10.4,"   max error = ",e10.5)', time, u_err/real(icount), u_maxerr
      end do
    end block
    close(fid)

    print *, 'analytical ', analytical_1d_temperature(width/2.0, t_end, &
        conductivity/capacity, width, temp_bc)

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
