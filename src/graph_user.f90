! -----------------------------------------------------------------------------
! User of graph library must set the storage space accordingly to application.
! These constants are used in "graph.f90".
!
! Also see "testsrc/map.f90" for mapping the storage space.
! -----------------------------------------------------------------------------
  module graph_user_mod
    implicit none (type, external)
    public

    ! Size of vertex/edge storage space
    integer, parameter :: VSIZE_IPAR = 4, VSIZE_RPAR = 7
    integer, parameter :: ESIZE_IPAR = 2, ESIZE_RPAR = 4
  end module graph_user_mod
