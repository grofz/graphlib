program test_mesh
  use iso_fortran_env, only : dp=>real64, output_unit
  use mesh_mod, only : mesh_t, mesh_handle_t
  use vtuio_mod, only : vtuio_write
  implicit none (type, external)

  type(mesh_t) :: grid, grid3d
  type(mesh_handle_t), allocatable :: parr(:), carr(:)
  integer :: i


  call grid%initialize(is_3d=.false.)

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

  call grid%print(output_unit)

  call vtuio_write('mesh', grid, [0,0,0,0])

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
  call grid3d%print(output_unit)
  call vtuio_write('mesh3d', grid3d, [0,0,0,0])

  print *, grid%npoints_per_cell(), 'expecting 3'
  print *, grid3d%npoints_per_cell(), 'expecting 4'
  print *, grid3d%graph_t%npoints_per_cell(), 'expecting 2'


end program
