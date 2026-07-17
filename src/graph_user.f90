  module graph_user_mod
    implicit none (type, external)
    public

    ! Size of vertex/edge storage space
    integer, parameter :: VSIZE_IPAR = 3, VSIZE_RPAR = 7
    integer, parameter :: ESIZE_IPAR = 1, ESIZE_RPAR = 2

    ! Maping of properties to the storeage space
    ! Vertex properties
    ! - integers
    integer, parameter :: VPOS_TYPE  = 1 ! kind of particle
    integer, parameter :: VPOS_BC    = 2 ! boundary condition
    integer, parameter :: VPOS_IWORK = 3 ! temporary working space
    ! - reals
    integer, parameter :: VPOS_RADIUS = 1 ! particle radius
    integer, parameter :: VPOS_X      = 2 ! position vector (3)
    integer, parameter :: VPOS_RWORK  = 5 ! temporary working space
    integer, parameter :: VPOS_RWORK2 = 6 ! temporary working space
    integer, parameter :: VPOS_RWORK3 = 7 ! temporary working space

    ! Edge properties
    ! - integers
    integer, parameter :: EPOS_TYPE = 1    ! edge type
    ! - reals
    integer, parameter :: EPOS_WEIGHT = 1  ! edge weight
    integer, parameter :: EPOS_RWORK = 2   ! temporary working space

    ! Mask for vtuio read/write
    integer, parameter :: &
        VTUIO_MASK(*) = [VPOS_RADIUS, VPOS_X, VPOS_TYPE, EPOS_TYPE]
  end module graph_user_mod
