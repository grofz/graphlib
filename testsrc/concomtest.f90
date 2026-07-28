program test_connected_components
  use graph_mod, only : graph_t
  use graph_user_mod
  use graph_testutils_mod, only : graph_from_arrays
  use vtuio_mod, only : vtuio_read, vtuio_write, vtuio_data_t
  use iso_fortran_env, only : output_unit
  implicit none (type, external)

  type(graph_t) :: g
  type(vtuio_data_t) :: vtudata
  integer :: ncomp
  integer, allocatable :: labels(:)

  ! ==================
  ! Test #1 a 2D plane
  ! ==================
  call vtuio_read('smallsample', g, VTUIO_MASK)
  
  g%vertices(1:g%nvertices)%ipar(VPOS_COMP) = -2 
  call g%connected_components( &
     !position_label=VPOS_COMP, &
      labels=labels, &
      vmask=g%vertices(1:g%nvertices)%ipar(VPOS_TYPE)==1, &
      lab_count=ncomp)
  g%vertices(1:g%nvertices)%ipar(VPOS_COMP)=labels
  print '("connented_components = ",i0)', ncomp
  block
    integer :: i
    do i=0, maxval(labels)
      print '("Value ",i0," is present ",i0," times in labels array")', i, count(labels==i)
    end do
  end block

  call vtudata%add_item('conlab', start=VPOS_COMP, iclass=0, ncomp=1, nbytes=4)
  call vtuio_write('a', g, VTUIO_MASK, vtudata=vtudata)

  ! ==========================
  ! Test #2 - undirected graph 
  ! ==========================
  200 block
    integer, allocatable :: cons(:,:)
    logical, allocatable :: vmask(:)
    integer, allocatable :: expected_labels(:)

    cons = reshape( [ &
      1, 2, &
      2, 7, &
      7, 1, &
      2, 3, &
      5, 3, &
      6, 4  &
      ], shape=[2,6])
    allocate(vmask(maxval(cons)), source=.true.)
    vmask(3) = .false.

    call graph_from_arrays(g, cons, is_directed_graph=.true.)
   !call graph_from_arrays(g, cons, is_directed_graph=.false.)

    expected_labels = [1, 1, 1, 2, 1, 2, 1]
    call g%connected_components(labels=labels)
    print '("Connected component labels ",*(i0,1x))', labels
    print '("Expecting ",*(i0,1x))', expected_labels
    print '("Passed = ",l1)', all(expected_labels==labels)

    expected_labels = [1, 1, 0, 2, 3, 2, 1]
    call g%connected_components(labels=labels, vmask=vmask)
    print '("Connected component labels ",*(i0,1x))', labels
    print '("Expecting ",*(i0,1x))', expected_labels
    print '("Passed = ",l1)', all(expected_labels==labels)
  end block
end program test_connected_components