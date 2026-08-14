! Copyright (C) 2026 Zdenek Grof
!
! This file is part of Graph library.
!
! Graph library is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published
! by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Graph library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with Graph library. If not, see <https://www.gnu.org/licenses/>.


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