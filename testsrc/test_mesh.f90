program test_mesh
  use iso_fortran_env, only : dp=>real64, output_unit
  use graph_mod, only : graph_handle_t => handle_t
  use mesh_mod, only : mesh_t
  use vtuio_mod, only : vtuio_write, vtuio_data_t, vtuio_read
  use map_mod, only : VPOS_TYPE, VPOS_X
  implicit none (type, external)

  type(vtuio_data_t) :: vtuio
  type(mesh_t) :: grid, grid3d, grid_read, grid3d_read
  type(graph_handle_t), allocatable :: parr(:), carr(:)
  integer :: i


  call grid%initialize(is_3d=.false.)
  call vtuio%add_item('type', VPOS_TYPE, 0, 1, 4)
  call vtuio%add_item('velo', VPOS_X, 2, 3, 8)
  call vtuio%add_item('ignore', 1, 3, 1, 4)

  allocate(parr(11))
  parr(1) =grid%add_point(real([1.0, 0.0, 0.0],dp))
  parr(2) =grid%add_point(real([3.0, 0.0, 0.0],dp))
  parr(3) =grid%add_point(real([5.0, 0.0, 0.0],dp))
  parr(4) =grid%add_point(real([7.0, 0.0, 0.0],dp))
  parr(5) =grid%add_point(real([2.0, 1.0, 0.0],dp))
  parr(6) =grid%add_point(real([4.0, 1.0, 0.0],dp))
  parr(7) =grid%add_point(real([6.0, 1.0, 0.0],dp))
  parr(8) =grid%add_point(real([1.0, 2.0, 0.0],dp))
  parr(9) =grid%add_point(real([3.0, 2.0, 0.0],dp))
  parr(10)=grid%add_point(real([5.0, 2.0, 0.0],dp))
  parr(11)=grid%add_point(real([7.0, 2.0, 0.0],dp))

  allocate(carr(10))
  carr(1)=grid%add_cell([parr(1),parr(2),parr(5), parr(1)])
  carr(2)=grid%add_cell([parr(2),parr(3),parr(6), parr(1)])
  carr(3)=grid%add_cell([parr(7),parr(3),parr(4), parr(1)])
  carr(4)=grid%add_cell([parr(1),parr(5),parr(8), parr(1)])
  carr(5)=grid%add_cell([parr(5),parr(6),parr(2), parr(1)])
  carr(6)=grid%add_cell([parr(5),parr(8),parr(9), parr(1)])
  carr(7)=grid%add_cell([parr(5),parr(6),parr(9), parr(1)])
  carr(8)=grid%add_cell([parr(6),parr(9),parr(10), parr(1)])
  carr(9)=grid%add_cell([parr(10),parr(11),parr(7), parr(1)])
  carr(10)=grid%add_cell([parr(4),parr(11),parr(7), parr(1)])

  do i=1, grid%nvertices
    grid%vertices(i)%ipar(VPOS_TYPE) = i
    grid%vertices(i)%rpar(VPOS_X:VPOS_X+2) = grid%points(i)%position
  end do
  call grid%print(output_unit)
  call vtuio_write('mesh', grid, vtudata=vtuio)

  ! Now test removing some cells
  call grid%remove_cell(carr(6))
  call grid%remove_cell(carr(5))
  call grid%remove_cell(carr(2))
  call grid%remove_cell(carr(7))
  call grid%remove_cell(carr(9))
  call grid%remove_cell(carr(1))
  call grid%remove_cell(carr(4))
  call grid%remove_cell(carr(8))
  call grid%remove_cell(carr(3))
  call grid%remove_cell(carr(10))
  print *, 'removed ok'
  call grid%print(output_unit)
  call vtuio_write('mesh1', grid, vtudata=vtuio)
  stop 1

  deallocate(parr, carr)
  call grid3d%initialize(is_3d=.true.)
  allocate(parr(6))
  parr(1)= grid3d%add_point([1.0_DP,0.0_DP,0.0_DP])
  parr(2)= grid3d%add_point([0.0_DP,1.0_DP,0.0_DP])
  parr(3)= grid3d%add_point([1.0_DP,1.0_DP,0.0_DP])
  parr(4)= grid3d%add_point([0.5_DP,0.5_DP,0.9_DP])
  parr(5)= grid3d%add_point([-0.5_DP,-0.5_DP,-0.5_DP])

  allocate(carr(4))
  carr(1) = grid3d%add_cell([parr(1),parr(2),parr(3),parr(4)])
  carr(1) = grid3d%add_cell([parr(1),parr(2),parr(3),parr(5)])
  carr(3) = grid3d%add_cell([parr(1),parr(2),parr(4),parr(5)])
  do i=1, grid%nvertices
    grid3d%vertices(i)%ipar(VPOS_TYPE) = i
  end do
  call grid3d%print(output_unit)
  call vtuio_write('mesh3d', grid3d, vtudata=vtuio)

  print *, 'TEST READ'
  call vtuio_read('mesh', grid_read, vtudata=vtuio)
  call grid_read%print(output_unit)
  call vtuio_write('mesh_copy', grid_read, vtudata=vtuio)

  call vtuio_read('mesh3d', grid3d_read, vtudata=vtuio)
  call grid3d_read%print(output_unit)
  call vtuio_write('mesh3d_copy', grid3d_read, vtudata=vtuio)
end program
