  program test_meshgen
    use graph_mod, only : graph_handle_t=>handle_t
    use mesh_mod, only : mesh_t, integrate_pde
    use iso_fortran_env, only : dp=>real64, output_unit
    implicit none (type, external)

    integer :: fid, i
    real(dp) :: conductivity, capacity
    real(dp), parameter :: h(*) = real([0.5,0.25,0.125,0.0625,0.03125],dp)
    real(dp), parameter :: dt(*)= real([1.0, 0.5, 0.25, 0.125, 0.0625],dp)

    interface
      subroutine make_timestep_study( &
          cell_size, dt_comps, conductivity, capacity, output_id)
        import dp
        implicit none (type, external)
        real(dp), intent(in) :: cell_size, dt_comps(:), conductivity, capacity
        integer, intent(in) :: output_id
      end subroutine
    end interface

    ! Parameters
    conductivity = 1.0_dp
    capacity = 1.0_dp
    open(newunit=fid, file='convergence11.log', status='replace')
    do i=1, size(h)
      call make_timestep_study( &
        cell_size = h(i), dt_comps = dt, &
        conductivity = 1.0_dp, capacity = 1.0_dp, output_id = fid)
      write(fid,*)
    end do
    close(fid)

    conductivity = 2.0_dp
    capacity = 1.0_dp
    open(newunit=fid, file='convergence21.log', status='replace')
    do i=1, size(h)
      call make_timestep_study( &
        cell_size = h(i), dt_comps = dt, &
        conductivity = 1.0_dp, capacity = 1.0_dp, output_id = fid)
      write(fid,*)
    end do
    close(fid)

    conductivity = 1.0_dp
    capacity = 2.0_dp
    open(newunit=fid, file='convergence12.log', status='replace')
    do i=1, size(h)
      call make_timestep_study( &
        cell_size = h(i), dt_comps = dt, &
        conductivity = 1.0_dp, capacity = 1.0_dp, output_id = fid)
      write(fid,*)
    end do
    close(fid)

!   call vtuio_write('meshgen', m, vtudata=vtudata)
!   call vtuio_write('meshgenv', m%graph_t, position_id=VPOS_X, vtudata=vtudata)
!   call vtuio_read('meshgen', m2, vtudata=vtudata)

  end program test_meshgen


  subroutine make_timestep_study(cell_size, dt_comps, conductivity, capacity, &
      output_id)
    use graph_mod, only : graph_handle_t=>handle_t
    use map_mod, only : VPOS_C, VPOS_X, VPOS_RTMP, VPOS_VB, EPOS_RTMP, VPOS_BC
    use mesh_mod, only : mesh_t, integrate_pde
    use vtuio_mod, only : vtuio_write, vtuio_data_t, vtuio_read
    use analytical_mod, only : analytical_1d_temperature
    use iso_fortran_env, only : dp=>real64
    implicit none (type, external)

    real(dp), intent(in) :: cell_size, dt_comps(:), conductivity, capacity
    integer, intent(in) :: output_id
!
! Test average error at two times for different time steps.
!
    type(mesh_t) :: m
    real(dp), parameter :: &
      width = 10.0_dp, &
      temp_bc(2) = [0.0_dp, 100.0_dp]
    real(dp) :: t_char, dt_comp, dt_out, t_end, t_start
    real(dp), allocatable :: u_init(:), u_out(:,:)
    integer, parameter :: NOUTS = 10
    integer, allocatable :: bc_label(:)
    integer :: irun

    t_char = capacity * width**2 / conductivity

    ! Generate mesh
    block
      real(dp), parameter :: rel_shift = 0.0_dp
      real(dp) :: p0(3), p1(3), p2(3)
      integer :: i, boffset(0:4)
      type(graph_handle_t), allocatable :: b(:)

      p0 = 0.0_dp
      p1 = [width, 0.0_dp, 0.0_dp]
      p2 = [0.0_dp, width, 0.0_dp]

      call m%initialize(is_3d=.false.)
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

      ! set conductivity
      m%vertices(1:m%nvertices)%rpar(VPOS_VB) = conductivity
      ! set capacity
      m%vertices(1:m%nvertices)%rpar(VPOS_RTMP) = capacity
    end block

    allocate(u_init(m%nvertices))

    900 format(2x,5(e12.5,1x))
    901 format(2x,5(a12,1x))
! 1234567890AB 1234567890AB 1234567890AB 1234567890AB 1234567890AB
    write(output_id,901) 'cell_size', 'dt_comp', 'time', 'u_err', 'u_maxerr'

    RUN_LOOP: do irun=1, size(dt_comps)
      ! Reset to the initial conditions
      where(bc_label==2)
        u_init = temp_bc(2)
      else where (bc_label==1)
        u_init = temp_bc(1)
      else where
        u_init = temp_bc(1)
      end where

      ! Integrate
      t_start = 0.0
      t_end = t_char
      dt_comp = dt_comps(irun)
      dt_out = t_char / real(NOUTS, dp)
      if (allocated(u_out)) deallocate(u_out)
      call integrate_pde(m, t_start, t_end, dt_comp, dt_out, &
        VPOS_VB, VPOS_RTMP, EPOS_RTMP, &
        bc_label, u_init, u_out)
      if (size(u_out,2)/=NOUTS+1) error stop 'size(u_out,2) unexpected'

      ! Compare numerical and analytical solution
      block
        integer :: u_count, itime(2), j, i
        real(dp) :: time, positions(3,4), y, u_anal, u_norm, u_err, u_maxerr

        itime = [2, size(u_out,2)]
        u_norm = (temp_bc(1) + temp_bc(2)) / 2.0_dp
        do j=1, 2
          time = t_start + min(t_end, real(itime(j)-1)*dt_out)
          u_err = 0.0
          u_maxerr = 0.0
          u_count = 0
          do i=1, m%nvertices
            if (m%v2c_index(i)<1) cycle
            positions = m%cells(m%v2c_index(i))%point_positions(m)
            y = sum(positions(2,1:3)) / 3.0_dp
           !y = m%vertices(i)%rpar(VPOS_X+1)
            u_anal = analytical_1d_temperature( &
                width-y, time, conductivity/capacity, width, temp_bc)
            u_err = u_err + abs(u_anal - u_out(i,itime(j)))
            u_maxerr = max(u_maxerr, abs(u_anal - u_out(i,itime(j))))
            u_count = u_count + 1
          end do
          u_err = u_err/real(u_count)
          write(output_id, 900) cell_size, dt_comp, time, u_err, u_maxerr
        end do
      end block

      ! write uout for the first tested t_out only
      if (irun==1) then
        block
          type(vtuio_data_t) :: vtudata
          integer :: i
          character(len=3) :: text

          ! What will be written to .vtu file for Paraview?
          call vtudata%add_item('conc', VPOS_C, 2, 1, 4)
          do i=1, size(u_out,2)
            write(text,'(i3.3)') i
            m%vertices(1:m%nvertices)%rpar(VPOS_C) = u_out(:,i)
            call vtuio_write('time'//text, m, vtudata=vtudata)
          end do
        end block
      end if
      flush(output_id)
    end do RUN_LOOP
  end subroutine make_timestep_study
