program test_mesh
  use iso_fortran_env, only : dp=>real64
  use mesh_mod, only : mesh_t, mesh_handle_t
  implicit none (type, external)

  type(mesh_t) :: grid
  type(mesh_handle_t), allocatable :: parr(:), carr(:)

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
end program