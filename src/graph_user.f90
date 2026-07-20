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
    integer, parameter :: VPOS_COMP  = 3 ! connected component
    ! - reals
    integer, parameter :: VPOS_RADIUS = 1 ! particle radius
    integer, parameter :: VPOS_X      = 2 ! position vector (2,3,4)
    integer, parameter :: VPOS_VB     = 5 ! vertex betweenness centrality
    integer, parameter :: VPOS_RWORK  = 6 ! temporary working space
    integer, parameter :: VPOS_C      = 7 ! concentration (potential)

    ! Edge properties
    ! - integers
    integer, parameter :: EPOS_TYPE = 1    ! edge type
    ! - reals
    integer, parameter :: EPOS_WEIGHT = 1  ! edge weight
    integer, parameter :: EPOS_EB     = 2  ! edge betweenness centrality
    integer, parameter :: EPOS_FLOW   = 3  ! flow

    ! Mask for vtuio read/write operation
    integer, parameter :: &
        VTUIO_MASK(*) = [VPOS_RADIUS, VPOS_X, VPOS_TYPE, EPOS_TYPE]
  end module graph_user_mod
