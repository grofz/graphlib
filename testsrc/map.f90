! -----------------------------------------------------------------------------
! Maping of properties to the storage space.
! Used by unit test programs.
!
! Also see "src/graph_user.f90" and verify that *SIZE_*PAR constants
! are large enough.
! -----------------------------------------------------------------------------
  module map_mod
    implicit none (type, external)
    public

    ! Vertex properties
    ! - integers (VSIZE_IPAR = 4)
    integer, parameter :: VPOS_TYPE  = 1 ! kind of particle
    integer, parameter :: VPOS_BC    = 2 ! boundary condition
    integer, parameter :: VPOS_COMP  = 3 ! connected component
    integer, parameter :: VPOS_ITMP  = 4 ! temporary working space
    ! - reals (VSIZE_RPAR = 7)
    integer, parameter :: VPOS_RADIUS = 1 ! particle radius
    integer, parameter :: VPOS_X      = 2 ! position vector (2,3,4)
    integer, parameter :: VPOS_VB     = 5 ! vertex betweenness centrality
    integer, parameter :: VPOS_RTMP   = 6 ! temporary working space
    integer, parameter :: VPOS_C      = 7 ! concentration (potential)

    ! Edge properties
    ! - integers (ESIZE_IPAR = 2)
    integer, parameter :: EPOS_TYPE = 1    ! edge type
    integer, parameter :: EPOS_ITMP = 2    ! temporary working space
    ! - reals (ESIZE_RPAR = 4)
    integer, parameter :: EPOS_WEIGHT = 1  ! edge weight
    integer, parameter :: EPOS_EB     = 2  ! edge betweenness centrality
    integer, parameter :: EPOS_FLOW   = 3  ! flow
    integer, parameter :: EPOS_RTMP   = 4  ! temporary working space

    ! Mask for vtuio read/write operation
    integer, parameter :: &
        VTUIO_MASK(*) = [VPOS_RADIUS, VPOS_X, VPOS_TYPE, EPOS_TYPE]
  end module map_mod
