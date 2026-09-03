! --------------------------------------------
! A unit test for cell_geometry_t calculations
! --------------------------------------------

  program cell_geometry_test
    use utest_mod, only : utest_t
    use graph_mod, only : graph_handle_t=>handle_t
    use mesh_mod, only : mesh_t, cell_geometry_t
    use vtuio_mod, only : vtuio_write
    use iso_fortran_env, only : dp=>real64
    implicit none(type, external)

    abstract interface
      type(cell_geometry_t) function exp_geometry_int(a)
        import dp, cell_geometry_t
        implicit none
        real(dp), intent(in) :: a
      end function
    end interface

    interface
      subroutine print_geometry(geom, expgeom)
        import cell_geometry_t
        implicit none
        type(cell_geometry_t), intent(in) :: geom, expgeom
      end subroutine
      subroutine build_geometry(positions, geometry_fun, a, geom, expgeom)
        import mesh_t, cell_geometry_t, graph_handle_t, dp, exp_geometry_int
        implicit none (type, external)
        real(dp), intent(in) :: positions(:,:)
        procedure(exp_geometry_int) :: geometry_fun
        real(dp), intent(in) :: a
        type(cell_geometry_t), intent(out) :: geom, expgeom
      end subroutine
      subroutine test_geometry(utest, geom, expgeom, msg, is_3d)
        import utest_t, cell_geometry_t, dp
        implicit none (type, external)
        type(utest_t), intent(inout) :: utest
        type(cell_geometry_t), intent(in) :: geom, expgeom
        character(len=*), intent(in) :: msg
        logical, intent(in) :: is_3d
      end subroutine
    end interface

    type(utest_t) :: utest
    type(mesh_t) :: mesh
    real(dp), allocatable :: positions(:,:)
    type(cell_geometry_t) :: geom, expgeom
    procedure(exp_geometry_int) :: expected_equilateral_triangle, &
        expected_right_triangle, expected_tetrahedron
    real(dp) :: rotation(3,3), alf
    real(dp), parameter :: a = 5.0_dp

    alf = -40.0_dp
    rotation = reshape( &
      [ cosd(alf), sind(alf), 0.0_dp, &
       -sind(alf), cosd(alf), 0.0_dp, &
       0.0_dp, 0.0_dp, 1.0_dp], shape=[3,3])

    ! equi-lateral triangle
    print '("EQUILATERAL TRIANGLE")'
    positions = reshape( &
      [ 0.0_dp, 0.0_dp, 0.0_dp, &
        a,      0.0_dp, 0.0_dp, &
        a/2,    sqrt(0.75)*a, 0.0_dp], shape=[3,3])
    call build_geometry(positions, expected_equilateral_triangle, a, &
        geom, expgeom)
    call print_geometry(geom, expgeom)
    call test_geometry(utest, geom, expgeom, 'eq-lateral triangle', .false.)

    ! rotate it
    print '("EQUILATERAL TRIANGLE ROTATED")'
    positions = matmul(rotation, positions)
    call build_geometry(positions, expected_equilateral_triangle, a, &
        geom, expgeom)
    expgeom%area_vector = matmul(rotation, expgeom%area_vector)
    expgeom%centre = matmul(rotation, expgeom%centre)
    call print_geometry(geom, expgeom)
    call test_geometry(utest, geom, expgeom, 'eq-lateral triangle rotated', .false.)

    ! right-angle triangle
    print '("RIGHT TRIANGLE")'
    positions = reshape( &
      [ 0.0_dp, 0.0_dp, 0.0_dp, &
        a,      0.0_dp, 0.0_dp, &
        0.0_dp, a,      0.0_dp], shape=[3,3])
    call build_geometry(positions, expected_right_triangle, a, &
        geom, expgeom)
    call print_geometry(geom, expgeom)
    call test_geometry(utest, geom, expgeom, 'right triangle', .false.)

    ! rotate it
    print '("RIGHT TRIANGLE ROTATED")'
    positions = matmul(rotation, positions)
    call build_geometry(positions, expected_right_triangle, a, &
        geom, expgeom)
    expgeom%area_vector = matmul(rotation, expgeom%area_vector)
    expgeom%centre = matmul(rotation, expgeom%centre)
    call print_geometry(geom, expgeom)
    call test_geometry(utest, geom, expgeom, 'right triangle rotated', .false.)

    ! regular tetrahedron
    print '("TETRAHEDRON")'
    positions = reshape( &
      [ 1.0_dp,  0.0_dp, -1.0_dp/sqrt(2.0_dp), &
       -1.0_dp,  0.0_dp, -1.0_dp/sqrt(2.0_dp), &
        0.0_dp,  1.0_dp,  1.0_dp/sqrt(2.0_dp), &
        0.0_dp, -1.0_dp,  1.0_dp/sqrt(2.0_dp) ], shape=[3,4])
    call build_geometry(positions, expected_tetrahedron, 2.0_dp, &
        geom, expgeom)
    call print_geometry(geom, expgeom)
    call test_geometry(utest, geom, expgeom, 'regular tetrahedron', .true.)

    ! rotate it
    print '("TETRAHEDRA ROTATED")'
    positions = matmul(rotation, positions)
    call build_geometry(positions, expected_tetrahedron, 2.0_dp, &
        geom, expgeom)
    expgeom%area_vector = matmul(rotation, expgeom%area_vector)
    expgeom%centre = matmul(rotation, expgeom%centre)
    call print_geometry(geom, expgeom)
    call test_geometry(utest, geom, expgeom, 'tetrahedron rotated', .true.)

   !call vtuio_write('cell_geometry', mesh)

    call utest%summarize()
    if (.not. utest%all_passed()) stop 1
  end program cell_geometry_test


  subroutine test_geometry(utest, geom, expgeom, msg, is_3d)
    use utest_mod, only : utest_t
    use mesh_mod, only : cell_geometry_t
    use iso_fortran_env, only : dp=>real64
    implicit none (type, external)
    type(utest_t), intent(inout) :: utest
    type(cell_geometry_t), intent(in) :: geom, expgeom
    character(len=*), intent(in) :: msg
    logical, intent(in) :: is_3d

    real(dp), parameter :: tol = 1.0e-7_dp
    integer :: n

    n = 3
    if (is_3d) n = 4

    call utest%within_tolerance(geom%volume, expgeom%volume, tol, &
      msg//': volume')
    call utest%within_tolerance(geom%area_vector(:,1), &
        expgeom%area_vector(:,1), tol, msg//': area vector 1')
    call utest%within_tolerance(geom%area_vector(:,2), &
        expgeom%area_vector(:,2), tol, msg//': area vector 2')
    call utest%within_tolerance(geom%area_vector(:,3), &
        expgeom%area_vector(:,3), tol, msg//': area vector 3')
    if (is_3d) &
      call utest%within_tolerance(geom%area_vector(:,4), &
          expgeom%area_vector(:,4), tol, msg//': area vector 4')
    call utest%within_tolerance(geom%face_distance(1:n), &
        expgeom%face_distance(1:n), tol, msg//': face distances')
    call utest%within_tolerance(geom%centre, expgeom%centre, tol, &
        msg//': centre')
  end subroutine test_geometry


  subroutine build_geometry(positions, exp_geometry_fun, a, geom, expgeom)
    use mesh_mod, only : mesh_t, cell_geometry_t
    use graph_mod, only : graph_handle_t=>handle_t
    use iso_fortran_env, only : dp=>real64
    implicit none (type, external)

    abstract interface
      type(cell_geometry_t) function exp_geometry_int(a)
        import dp, cell_geometry_t
        implicit none
        real(dp), intent(in) :: a
      end function
    end interface

    real(dp), intent(in) :: positions(:,:)
    procedure(exp_geometry_int) :: exp_geometry_fun
    real(dp), intent(in) :: a
    type(cell_geometry_t), intent(out) :: geom, expgeom

    type(mesh_t) :: mesh
    type(graph_handle_t), allocatable :: points(:), cells(:)
    type(graph_handle_t) :: null_handle
    integer :: i

    if (size(positions,2)==3) then
      call mesh%initialize(is_3d = .false.)
    else if (size(positions,2)==4) then
      call mesh%initialize(is_3d = .true.)
    else
      error stop 'build_geometry - 3 or 4 points expected'
    end if
    allocate(points(size(positions,2)))
    do i=1, size(points)
      points(i) = mesh%add_point(positions(:,i))
    end do

    if (allocated(cells)) deallocate(cells)
    allocate(cells(1))
    if (mesh%is_3d()) then
      cells(1) = mesh%add_cell(points(1:4))
    else
      cells(1) = mesh%add_cell([points(1:3), null_handle])
    end if
    geom = mesh%cells(mesh%index_from_handle(cells(1)))%geometry(mesh)
    expgeom = exp_geometry_fun(a)
  end subroutine build_geometry


  subroutine print_geometry(geom, expgeom)
    use mesh_mod, only : cell_geometry_t
    implicit none (type, external)
    type(cell_geometry_t), intent(in) :: geom, expgeom

    integer :: i

    111 format(SP, (a10, 2x, 4(g0.6,3x)))
    print 111, 'cent', geom%centre
    print 111, 'cent exp', expgeom%centre
    print 111, 'volu', geom%volume
    print 111, 'volu exp', expgeom%volume
    print 111, ( &
      'avec', geom%area_vector(:,i), norm2(geom%area_vector(:,i)), &
      'avec exp', expgeom%area_vector(:,i), norm2(expgeom%area_vector(:,i)), &
      i=1,4)
    print 111, 'fdis', (geom%face_distance(i), i=1,4)
    print 111, 'fdis exp', (expgeom%face_distance(i), i=1,4)
  end subroutine print_geometry


  type(cell_geometry_t) function expected_equilateral_triangle(a) result(geom)
    use mesh_mod, only : cell_geometry_t
    use iso_fortran_env, only : dp => real64
    implicit none (type, external)
    real(dp), intent(in) :: a

    geom%centre = [a/2.0_dp, sqrt(0.75_dp)/3.0_dp*a, 0.0_dp]
    geom%volume = 0.5_dp*sqrt(0.75_dp)*a**2
    geom%face_distance = sqrt(0.75_dp)/3.0_dp*a       ! same for all 3-faces
    geom%area_vector(:,1) = [+a*cosd(30.0_dp), a*sind(30.0_dp), 0.0_dp]
    geom%area_vector(:,2) = [-a*cosd(30.0_dp), a*sind(30.0_dp), 0.0_dp]
    geom%area_vector(:,3) = [0.0_dp,        -a,           0.0_dp]
  end function


  type(cell_geometry_t) function expected_right_triangle(a) result(geom)
    use mesh_mod, only : cell_geometry_t
    use iso_fortran_env, only : dp => real64
    implicit none (type, external)
    real(dp), intent(in) :: a

    geom%centre = [a/2.0_dp, a/2.0_dp, 0.0_dp]
    geom%volume = 0.5_dp*a**2
    geom%face_distance(1) = 0
    geom%face_distance(2) = 0.5_dp*a
    geom%face_distance(3) = 0.5_dp*a

    geom%area_vector(:,1) = [a, a, 0.0_dp]
    geom%area_vector(:,2) = [-a, 0.0_dp, 0.0_dp]
    geom%area_vector(:,3) = [0.0_dp, -a, 0.0_dp]
  end function


  type(cell_geometry_t) function expected_tetrahedron(a) result(geom)
    use mesh_mod, only : cell_geometry_t
    use iso_fortran_env, only : dp => real64
    implicit none (type, external)
    real(dp), intent(in) :: a

    ! Circumcentre is in [0, 0, 0]
    geom%centre = 0.0_dp

    ! Volume of tetrahedron with edge length "a" is
    geom%volume = a**3 / (6.0_dp*sqrt(2.0_dp))

    ! The area of each side is S = sqrt(3)/4 * a**2
    geom%area_vector(:,1) = [ -sqrt(2.0_dp), 0.0_dp, 1.0_dp]
    geom%area_vector(:,2) = [ 0.0_dp, -sqrt(2.0_dp), -1.0_dp]
    geom%area_vector(:,3) = [ sqrt(2.0_dp), 0.0_dp, 1.0_dp]
    geom%area_vector(:,4) = [ 0.0_dp, sqrt(2.0_dp), -1.0_dp]
    geom%area_vector = a**2/4.0_dp * geom%area_vector

    ! The radii of insphere (https://en.wikipedia.org/wiki/Regular_tetrahedron)
    ! is r_in = 1/3 * r_ci
    ! and r_ci = a * sqrt(3) / (2*sqrt(2))
    ! so r_in = a / (2*sqrt(6))
    geom%face_distance = a / (2.0_dp*sqrt(6.0_dp))
  end function
