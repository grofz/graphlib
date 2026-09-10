  program test_meshgen
    use iso_fortran_env, only : dp=>real64, output_unit
    implicit none (type, external)

    integer :: fid
    real(dp) :: conductivity, capacity
   !real(dp), parameter :: h(*) = real([0.5,0.25,0.125,0.0625,0.03125],dp)
    real(dp), parameter :: h(*) = real([0.5,0.25,0.125,0.0625],dp)
    real(dp), parameter :: dt(*)= real([1.0, 0.5, 0.25, 0.125, 0.0625],dp)

    interface
      subroutine make_timestep_study( &
          cell_sizes, dt_comps, conductivity, capacity, output_id)
        import dp
        implicit none (type, external)
        real(dp), intent(in) :: cell_sizes(:), dt_comps(:), conductivity, capacity
        integer, intent(in) :: output_id
      end subroutine
    end interface

    ! -
    conductivity = 1.0_dp
    capacity = 1.0_dp
    open(newunit=fid, file='convergence11.log', status='replace')
    call make_timestep_study( &
      cell_sizes = h, dt_comps = dt, &
      conductivity = conductivity, capacity = capacity, output_id = fid)
    close(fid)

    conductivity = 2.0_dp
    capacity = 1.0_dp
    open(newunit=fid, file='convergence21.log', status='replace')
    call make_timestep_study( &
      cell_sizes = h, dt_comps = dt, &
      conductivity = conductivity, capacity = capacity, output_id = fid)
    close(fid)

    conductivity = 1.0_dp
    capacity = 2.0_dp
    open(newunit=fid, file='convergence12.log', status='replace')
    call make_timestep_study( &
      cell_sizes = h, dt_comps = dt, &
      conductivity = conductivity, capacity = capacity, output_id = fid)
    close(fid)

  end program test_meshgen


  subroutine make_timestep_study(cell_sizes, dt_comps, conductivity, capacity, &
      output_id)
    use graph_mod, only : graph_handle_t=>handle_t
    use map_mod, only : VPOS_C, VPOS_X, VPOS_RTMP, VPOS_VB, EPOS_RTMP
    use mesh_mod, only : mesh_t, integrate_pde
    use vtuio_mod, only : vtuio_write, vtuio_data_t
    use analytical_mod, only : analytical_1d_temperature
    use iso_fortran_env, only : dp=>real64
    implicit none (type, external)

    real(dp), intent(in) :: cell_sizes(:), dt_comps(:), conductivity, capacity
    integer, intent(in) :: output_id
!
! Test average error at two times for different spatial and time steps.
!
! "collected_err" array maping:
!   1. index = cell size
!   2. index = time step
!   3. index = at 10% of char. time, at 100% of char. time
!   4. index = avg/max error
!
    real(dp) :: collected_err(size(cell_sizes), size(dt_comps), 2, 2)
    real(dp) :: t_char, dt_out, t_end, t_start, times(2)
    real(dp), allocatable :: u_init(:), u_out(:,:), t_out(:), cell_y(:)
    real(dp), parameter :: &
      width = 10.0_dp, &
      u_bc(2) = [0.0_dp, 100.0_dp]
    integer :: imesh, idt, i, itime(2)
    integer, allocatable :: bc_label(:)
    integer, parameter :: NOUTS = 10, ID_AVGERR=1, ID_MAXERR=2
    type(mesh_t) :: m

    write(output_id,'("capacity = ",g0,"   conductivity = ",g0)') &
        capacity, conductivity
    t_char = capacity * width**2 / conductivity
    t_start = 0.0
    t_end = t_start + t_char
    dt_out = t_char / real(NOUTS, dp)
    itime = [1+1, NOUTS+1]
    times(1) = t_start + 1.0_dp/real(NOUTS)*t_char
    times(2) = t_start + t_char

    MESH_LOOP: do imesh=1, size(cell_sizes)

      ! Generate mesh
      block
        real(dp), parameter :: rel_shift = 0.0_dp
        real(dp) :: p0(3), p1(3), p2(3), positions(3,4)
        integer :: boffset(0:4)
        type(graph_handle_t), allocatable :: b(:)

        p0 = 0.0_dp
        p1 = [width, 0.0_dp, 0.0_dp]
        p2 = [0.0_dp, width, 0.0_dp]

        call m%initialize(is_3d=.false.)
        call m%append_rectilinear_mesh( &
            p0, p1, p2, cell_sizes(imesh), rel_shift, VPOS_X, b, boffset)

        ! remove 02 and 13 ghost nodes (west / east)
        do i=boffset(1)+1, boffset(2)
          call m%remove_vertex(b(i))
        end do
        do i=boffset(3)+1, boffset(4)
          call m%remove_vertex(b(i))
        end do

        ! set boundary and initial conditions
        if (allocated(bc_label)) deallocate(bc_label)
        allocate(bc_label(m%nvertices), source=0)
        do i=boffset(2)+1, boffset(3) ! north border
          bc_label(m%index_from_handle(b(i))) = 2 ! u_2
        end do
        do i=boffset(0)+1, boffset(1) ! south border
          bc_label(m%index_from_handle(b(i))) = 1 ! u_1
        end do

        ! set conductivity
        m%vertices(1:m%nvertices)%rpar(VPOS_VB) = conductivity
        ! set capacity
        m%vertices(1:m%nvertices)%rpar(VPOS_RTMP) = capacity

        ! cell y-coordinate
        call m%build_v2c_index()
        if (allocated(cell_y)) deallocate(cell_y)
        allocate(cell_y(m%nvertices))
        do i=1, m%nvertices
          if (m%v2c_index(i)<1) cycle
          positions = m%cells(m%v2c_index(i))%point_positions(m)
          cell_y(i) = sum(positions(2,1:3)) / 3.0_dp
        end do

        print '("Mesh with ",i0," cells and ",i0," vertices generated")', &
            m%ncells, m%nvertices
      end block

      if (allocated(u_init)) deallocate(u_init)
      allocate(u_init(m%nvertices))


      RUN_LOOP: do idt=1, size(dt_comps)
        ! Reset to the initial conditions
        where(bc_label==2)
          u_init = u_bc(2)
        else where (bc_label==1)
          u_init = u_bc(1)
        else where
          u_init = u_bc(1)
        end where

        ! Integrate
        if (allocated(u_out)) deallocate(u_out)
        if (allocated(t_out)) deallocate(t_out)
        call integrate_pde(m, t_start, t_end, dt_comps(idt), dt_out, &
          VPOS_VB, VPOS_RTMP, EPOS_RTMP, &
          bc_label, u_init, t_out, u_out)
        if (size(u_out,2)/=NOUTS+1) error stop 'size(u_out,2) unexpected'

        ! Compare numerical and analytical solution
        block
          integer :: u_count, j, i
          real(dp) :: u_anal, u_err, u_maxerr

          do j=1, 2
            u_err = 0.0
            u_maxerr = 0.0
            u_count = 0
            do i=1, m%nvertices
              if (m%v2c_index(i)<1) cycle
              u_anal = analytical_1d_temperature( &
                  width-cell_y(i), times(j), &
                  conductivity/capacity, width, u_bc)
              u_err = u_err + abs(u_anal - u_out(i,itime(j)))
              u_maxerr = max(u_maxerr, abs(u_anal - u_out(i,itime(j))))
              u_count = u_count + 1
            end do
            u_err = u_err/real(u_count)
            collected_err(imesh, idt, j, ID_AVGERR) = u_err
            collected_err(imesh, idt, j, ID_MAXERR) = u_maxerr
          end do
          flush(output_id)
        end block

        ! write uout for the last "dt" and last "h" value only
        if (idt==size(dt_comps) .and. imesh==size(cell_sizes)) then
          block
            type(vtuio_data_t) :: vtudata
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
      end do RUN_LOOP
    end do MESH_LOOP

    ! print results in more convenient order than calculation order
    do i = 1, 2
      ! i=1: transient (at 10% of char. time)
      ! i=2: almost steady state (at 100% of char. time)
      select case(i)
      case(1)
        write(output_id,'("transient error effect")')
      case(2)
        write(output_id,'("steady-state error effect")')
      end select
      900 format(2x,5(e12.5,1x))
      901 format(2x,5(a12,1x))
      write(output_id,901) 'cell_size', 'dt_comp', 'time', 'u_avgerr', 'u_maxerr'
      do imesh=1, size(cell_sizes)
        do idt=1, size(dt_comps)
          write(output_id, 900) cell_sizes(imesh), dt_comps(idt), times(i), &
            collected_err(imesh, idt, i, ID_AVGERR), &
            collected_err(imesh, idt, i, ID_MAXERR)
        end do
        write(output_id,*)
      end do
    end do
  end subroutine make_timestep_study
