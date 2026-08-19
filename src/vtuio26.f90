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


! CONTENT (only public procedures are listed)
! module vtuio_mod
!   subroutine vtuio_write(file, graph, mask, time, vtudata)
!   subroutine vtuio_read(file, graph, mask, time)
!
! THIS IS A MODIFIED VERSION (July 2026)

  module vtuio_mod
    use graph_mod, only : graph_t, handle_t, edge_t, vertex_t
    use graph_user_mod, only : VSIZE_RPAR, VSIZE_IPAR, ESIZE_RPAR, ESIZE_IPAR
    use mesh_mod, only : mesh_t, point_t, cell_t
    use vtuio_tree_mod, only : object_t, vtuio_tree_read
    use iso_fortran_env, only : &
    &   SP=>real32, DP=>real64, I4B=>int32, I8B=>int64, I1B=>int8, &
    &   error_unit
    implicit none
    private
    public vtuio_write, vtuio_read

    type :: vtuio_meta_t
      integer :: start
        !! row index where to start in data-array
      integer :: ncomp
        !! number of components (1=scalar, 3=vector)
      integer :: nbytes
        !! number of bytes stored (1,4 for integers / 4,8 for reals)
      character(len=20) :: label
        !! label as shown in paraview
      integer(I1B) :: iclass
        !* -----------------------------
        !  | iclass   | points | cells |
        !  | -------- | ---------------|
        !  | integers |    0   |   1   |
        !  | reals    |    2   |   3   |
        !  -----------------------------
    end type vtuio_meta_t

    integer(I1B), parameter, public :: &
        META_IS_INT   = 0, META_IS_REAL  = 2, & ! 00_binary or 10_binary
        META_IS_POINT = 0, META_IS_CELL  = 1    ! 00_binary or 01_binary
      ! sum of INT/REAL and POINT/CELL options gives 0,1,2,3

    type, public :: vtuio_data_t
      !! 1. Call `add_item` to add additional data fields
      !! 2. Call vtuio_write
      !! 3. Call `finalize`
      type(vtuio_meta_t), allocatable :: meta(:)
    contains
      procedure :: add_item => meta_add_item
      procedure :: free => meta_free
    end type vtuio_data_t

    interface reallocate
      module procedure reallocate_real, reallocate_int
    end interface

    interface parse_value
      module procedure parse_value1, parse_value2
    end interface

! temporary to keep old interface available
interface vtuio_write
  module procedure vtuio_write1, vtuio_write2
end interface
interface vtuio_read
  module procedure vtuio_read1, vtuio_read2
end interface

! -----------------------------------------------------------------------------
! TUNING THE BINARY FORMAT FOR vtuio_write
!
! In AppendedData section, every DataArray piece starts with an integer
! stating the size of that piece in bytes. Integer used in binary data
! headers can be 32bit or 64bit (4=UInt32, 8=UInt64)
!
! 64bit option:
    integer, parameter :: HEADERTYPE_KIND=I8B, HEADERTYPE_SIZE=8
    character(len=*), parameter :: HEADERTYPE_TEXT='header_type="UInt64"'
! 32 bit option:
    !   integer, parameter :: header_type_kind=I4B, header_type_size=4
    !   character(len=*), parameter :: header_type_text='header_type="UInt32"'
!
! POINT POSITIONS:
! Float32
    integer, parameter :: POSITIONS_KIND=SP, POSITIONS_SIZE=4
    character(len=*), parameter :: POSITIONS_TEXT='"Float32"'
! Float64
    !   integer, parameter :: positions_kind=dp, positions_size=8
    !   character(len=*), parameter :: positions_text='"Float64"'
!
! CONNECTIVITY:
! Int32
    integer, parameter :: CONNECTIONS_KIND=I4B, CONNECTIONS_SIZE=4
    character(len=*), parameter :: CONNECTIONS_TEXT='"Int32"'
! Int64
    !   integer, parameter :: connections_kind=8, connections_size=8
    !   character(len=*), parameter :: connections_text='"Int64"'
    !
! TYPE OF CELLS (all cells have the same type)
    integer, parameter :: CELLTYPE_KIND=I1B, CELLTYPE_SIZE=1
    character(len=*), parameter :: CELLTYPE_TEXT='"UInt8"'
    integer, parameter :: VTK_LINE = 3, VTK_TRIANGLE = 5, VTK_TETRA = 10
! -----------------------------------------------------------------------------

    integer, parameter :: DATAARRAY_NOT_FOUND=-1
    integer, parameter :: MAX_BUFFER_LEN=400
    character(len=1), parameter :: LF=char(10) ! end-of-line
    character(len=*), parameter :: SUFFIX = '.vtu'

!TODO to be removedd
    ! legend for items in "mask" array argument
integer, parameter :: MASK_RADIUS=1, MASK_POSITION=2, MASK_POINT_TYPE=3, MASK_CELL_TYPE=4

  contains

    ! ----------------------
    ! Helper data descriptor
    ! ----------------------

    pure subroutine data_descriptors(int_or_real, nbytes, data_kind, data_text)
      integer, intent(in) :: int_or_real  ! IS_INT or IS_REAL
      integer, intent(in) :: nbytes       ! 1/4/8 for integers, 4/8 for reals
      integer, intent(out) :: data_kind
      character(len=10), intent(out) :: data_text

      select case(int_or_real)
      case(META_IS_REAL)
        select case(nbytes)
        case(4)
          data_kind = SP
          data_text = '"Float32"'
        case(8)
          data_kind = DP
          data_text = '"Float64"'
        case default
          error stop 'data_descriptors - 4 or 8 bytes for real data'
        end select
      case(META_IS_INT)
        select case(nbytes)
        case(1)
          data_kind = I1B
          data_text = '"UInt8"' ! unsigned integer!
        case(4)
          data_kind = I4B
          data_text = '"Int32"'
        case(8)
          ! Note: as our int arrays are only 4B, this is not necessary
          data_kind = I8B
          data_text = '"Int64"'
        case default
          error stop 'write_data - 1 4 8 bytes for int data'
        end select
      end select
    end subroutine data_descriptors


    pure subroutine get_counts(graph, npoints, ncells)
      class(graph_t), intent(in) :: graph
      integer, intent(out) :: npoints, ncells
!
! Mesh - points are points
!      - cells are cells (triangles or tetrahedra)      
!      - vertex data become cell data
!      - edge data is ignored
!
! Graph - vertices are points
!       - edges are cells (lines)
!       - vertex data become point data
!       - edge data become cell data
!
      select type(graph)
      class is (mesh_t)
        npoints = graph%npoints
        ncells = graph%ncells
        if (graph%ncells /= graph%nvertices) error stop &
            'vtuio_write - npoints == nvertices required'
        ! TODO (implementation progress note)
        ! Required at this moment.
        ! In future version also the presence of "non-mesh" vertices
        ! will be supported (e.g. ghost cells). 
        ! The requirement will be graph%ncells <= graph%nvertices
      class default
        npoints = graph%nvertices
        ncells = graph%nedges
      end select
    end subroutine get_counts


    ! -----------------------------
    ! Writing the Unstructured Grid
    ! -----------------------------
! temporary old interface
subroutine vtuio_write2(file, graph, mask, time, vtudata)
  character(len=*), intent(in) :: file
  class(graph_t), intent(in) :: graph
  integer, intent(in) :: mask(4)
    !! mask(1) - index pointing to "radius" in "vertex%rpar" array
    !! mask(2) - index pointing to the first "position" component in "vertex%rpar" array
    !! mask(3) - index pointing to "type" in "vertex%ipar" array
    !! mask(4) - index pointing to "type" in "edge%ipar" array
  real(DP), intent(in), optional :: time
  type(vtuio_data_t), optional :: vtudata
  type(vtuio_data_t) :: vtudata0

  print *, 'WARNING - using obsolete vtuio_write interface'
  if (present(vtudata)) then
    call vtudata%add_item('radius', mask(MASK_RADIUS), int(META_IS_REAL+META_IS_POINT), 1, 4)
    call vtudata%add_item('type', mask(MASK_POINT_TYPE), int(META_IS_INT+META_IS_POINT), 1, 1)
    call vtudata%add_item('con t', mask(MASK_CELL_TYPE), int(META_IS_INT+META_IS_CELL), 1, 1)
    call vtuio_write1(file, graph, mask(MASK_POSITION), time, vtudata)
  else
    call vtudata0%add_item('radius', mask(MASK_RADIUS), int(META_IS_REAL+META_IS_POINT), 1, 4)
    call vtudata0%add_item('type', mask(MASK_POINT_TYPE), int(META_IS_INT+META_IS_POINT), 1, 1)
    call vtudata0%add_item('con t', mask(MASK_CELL_TYPE), int(META_IS_INT+META_IS_CELL), 1, 1)
    call vtuio_write1(file, graph, mask(MASK_POSITION), time, vtudata0)
  end if
end subroutine

!TODO to be renamed to vtuio_write
    subroutine vtuio_write1(file, graph, position_id, time, vtudata)
      !* Write Unstructured Grid - vertices and edges
      character(len=*), intent(in) :: file
        !! file name without .vtu suffix
      class(graph_t), intent(in) :: graph
        !! contains vertices/points and edges/cells
      integer, intent(in), optional :: position_id
      real(DP), intent(in), optional :: time
      type(vtuio_data_t), optional :: vtudata
        !! data structure for additional data
!
! TODO Documentation Block
!
      character(len=MAX_BUFFER_LEN/2) :: text1, text2
      integer :: npoints, ncells, fid, offset

      if (.not. graph%is_initialized()) error stop &
          'vtuio_write - graph is not initialized'

      call get_counts(graph, npoints, ncells)
      offset = 0
      open(newunit=fid, file=trim(file)//SUFFIX, status='replace', &
      &    access='stream')

      write(fid) '<VTKFile type="UnstructuredGrid" version="1.0"&
      & byte_order="LittleEndian" ', HEADERTYPE_TEXT, '>', LF
      write(fid) '<UnstructuredGrid>', LF

      if (present(time)) call write_time_value(fid, time)

      write(text1,'(i0)') npoints
      write(text2,'(i0)') ncells
      write(fid) '<Piece NumberOfPoints="', trim(text1), &
      &          '" NumberOfCells="', trim(text2), '">', LF

      ! Point data - header
      write(fid) '  <PointData>', LF
      call extract_and_write_data(META_IS_POINT, &
          graph, fid, vtudata=vtudata, offset=offset)
      write(fid) '  </PointData>', LF

      ! Cell data - header
      write(fid) '  <CellData>', LF
      call extract_and_write_data(META_IS_CELL, &
          graph, fid, vtudata=vtudata, offset=offset)
      write(fid) '  </CellData>', LF

      ! Points - header
      write(fid) '  <Points>', LF
      call write_points(fid, graph, npoints, position_id, offset)
      write(fid) '  </Points>', LF

      ! Cells - header
      write(fid) '  <Cells>', LF
      call write_connectivity(fid, graph, ncells, offset)
      write(fid) '  </Cells>', LF

      ! Footer / Appended data
      write(fid) '</Piece>', LF
      write(fid) '</UnstructuredGrid>', LF
      write(fid) '<AppendedData encoding="raw">', LF
      write(fid) '_' ! data block starts with underscore

      ! Point data - binary
      call extract_and_write_data(META_IS_POINT, graph, fid, vtudata=vtudata)

      ! Cell data - binary
      call extract_and_write_data(META_IS_CELL, graph, fid, vtudata=vtudata)

      ! Points and cells - binary
      call write_points(fid, graph, npoints, position_id)
      call write_connectivity(fid, graph, ncells)

      ! Closing tags
      write(fid) ' ', LF
      write(fid) '</AppendedData>', LF
      write(fid) '</VTKFile>'
      close(fid)

      print '("VTU write: ",i0," points and ",i0," cells written&
          & to """,a,"""")', npoints, ncells, trim(file)//SUFFIX
    end subroutine vtuio_write1


    subroutine write_points(fid, graph, npoints, positions_id, offset)
      integer, intent(in) :: fid, npoints
      integer, intent(in), optional ::  positions_id
      class(graph_t), intent(in) :: graph
      integer, intent(inout), optional :: offset
!
! Called twice:
!   1. with offset present: write XML header and update offset
!   2. witout offset: write binary data 
! Calls must be co-ordinated with calls to write_data and write_connectivity.
!
      integer :: i
      character(len=MAX_BUFFER_LEN/2) :: text1

      if (present(offset)) then
        write(text1,'(i0)') offset
        write(fid) '    <DataArray type=', POSITIONS_TEXT, &
        & ' NumberOfComponents="3" format="appended" offset="', trim(text1), &
        & '" />', LF
        offset = offset + HEADERTYPE_SIZE + 3*npoints*POSITIONS_SIZE

      else
        write(fid) int(3*npoints*POSITIONS_SIZE, kind=HEADERTYPE_KIND)

        select type (graph)
        class is (mesh_t)
          do i=1, npoints
            write(fid) real(graph%points(i)%position, kind=POSITIONS_KIND)
          end do
        class default
          if (.not. present(positions_id)) error stop &
              'vtuio_write - position_id is a compulsory argument for graph_t'
          if (positions_id < 1 .or. positions_id+2 > VSIZE_RPAR) error stop &
              'vtuio_write - position_id is invalid'
          do i=1, npoints
            write(fid) real(graph%vertices(i)%rpar(positions_id:positions_id+2), &
                kind=POSITIONS_KIND)
          end do
        end select
      end if
    end subroutine write_points


    subroutine write_connectivity(fid, graph, ncells, offset)
      integer, intent(in) :: fid, ncells
      class(graph_t), intent(in) :: graph
      integer, intent(inout), optional :: offset
!
! Called twice:
!   1. with offset present: write XML header and update offset
!   2. witout offset: write binary data 
! The calls must be co-ordinated with calls to write_data and write_points.
!
      integer :: i, n, pids(4), id_ubound
      character(len=MAX_BUFFER_LEN/2) :: text1

      n = graph%npoints_per_cell()

      if (present(offset)) then
        ! Write XML header
        ! connectivity
        write(text1,'(i0)') offset
        write(fid) '    <DataArray type=', CONNECTIONS_TEXT, &
        & ' Name="connectivity" format="appended" offset="', trim(text1), &
        & '" />', LF
        offset = offset + HEADERTYPE_SIZE + n*ncells*CONNECTIONS_SIZE

        ! offsets
        write(text1,'(i0)') offset
        write(fid) '    <DataArray type=', CONNECTIONS_TEXT, &
        & ' Name="offsets" format="appended" offset="', trim(text1), '" />', LF
        offset = offset + HEADERTYPE_SIZE + ncells*CONNECTIONS_SIZE

        ! types
        write(text1,'(i0)') offset
        write(fid) '    <DataArray type=', CELLTYPE_TEXT, &
        & ' Name="types" format="appended" offset="', trim(text1), '" />', LF
        offset = offset + HEADERTYPE_SIZE + ncells*CELLTYPE_SIZE

      else
        ! Write binary data
        select type(graph)
        class is (mesh_t)
          id_ubound = graph%npoints
        class default
          id_ubound = graph%nvertices
        end select

        ! connectivity
        write(fid) int(n*ncells*CONNECTIONS_SIZE, kind=HEADERTYPE_KIND)

        do i=1, ncells
          select type (graph)
          class is (mesh_t)
            pids = graph%cells(i)%point_indices(graph)
          class default
            pids(1:n) = graph%edges(i)%vertex_indices(graph)
          end select
          if (any(pids(1:n)<1) .or. any(pids(1:n)>id_ubound)) &
            &    error stop 'write_connectivity - cell references an unknown vertex/point'

          ! for Paraview, index starts at zero (not at one)
          select case(n)
          case(3) ! 2D-mesh (triangles)
            write(fid) &
                int(pids(1)-1, kind=CONNECTIONS_KIND), &
                int(pids(2)-1, kind=CONNECTIONS_KIND), &
                int(pids(3)-1, kind=CONNECTIONS_KIND)
          case(4) ! 3D-mesh (tetrahedra)
            write(fid) &
                int(pids(1)-1, kind=CONNECTIONS_KIND), &
                int(pids(2)-1, kind=CONNECTIONS_KIND), &
                int(pids(3)-1, kind=CONNECTIONS_KIND), &
                int(pids(4)-1, kind=CONNECTIONS_KIND)
          case(2) ! vertices/edges (lines)
            write(fid) &
                int(pids(1)-1, kind=CONNECTIONS_KIND), &
                int(pids(2)-1, kind=CONNECTIONS_KIND)
          case default
            error stop 'vtuio_write - npoints_per_cell invalid (internal error'
          end select
        end do

        ! offsets (all cells consist of the same number of points)
        write(fid) int(ncells*CONNECTIONS_SIZE, kind=HEADERTYPE_KIND)
        write(fid) (int(i*n, kind=CONNECTIONS_KIND), i=1,ncells)

        ! types (all cells are of the same type)
        write(fid) int(ncells*CELLTYPE_SIZE, kind=HEADERTYPE_KIND)
        select case (n)
        case(3)
          write(fid) (int(VTK_TRIANGLE, kind=CELLTYPE_KIND), i=1,ncells)
        case(4)
          write(fid) (int(VTK_TETRA, kind=CELLTYPE_KIND), i=1,ncells)
        case(2)
          write(fid) (int(VTK_LINE, kind=CELLTYPE_KIND), i=1,ncells)
        end select
      end if
    end subroutine write_connectivity


    subroutine write_data(fid, nbytes, label, rdata, idata, offset)
      integer, intent(in) :: fid
      integer, intent(in) :: nbytes ! 1, 4, 8
      character(len=*), intent(in) :: label
      real(DP), intent(in), optional :: rdata(:,:)
      integer, intent(in), optional :: idata(:,:)
      integer, intent(inout), optional :: offset
!
! Called twice for each DataArray:
!   1. with offset present: write XML header and update offset
!   2. witout offset: write binary data 
! The calls must occur in exactly same order and co-ordinated with calls
! to write_connectivity and write_points.
!
! Array "rdata" or "idata" must be present (not both)
! - dimension 1 is component (scalars or vectors)
! - dimension 2 is point / connection
! If "offset" is present then write the header, otherwise write data
!
      integer :: real_or_int ! IS_REAL or IS_INT
      integer :: data_kind, ncomps, n
      character(len=10) :: data_text
      character(len=MAX_BUFFER_LEN/2) :: text1, text2

      if (present(rdata) .and. .not. present(idata)) then
        real_or_int = META_IS_REAL
      else if (.not. present(rdata) .and. present(idata)) then
        real_or_int = META_IS_INT
      else
        error stop 'write_data - rdata/idata must be present, but not both'
      end if

      call data_descriptors(real_or_int, nbytes, data_kind, data_text)

      select case(real_or_int)
      case(META_IS_REAL)
        ncomps = size(rdata, dim=1)
        n = size(rdata, dim=2)
      case(META_IS_INT)
        ncomps = size(idata, dim=1)
        n = size(idata, dim=2)
      case default
        error stop 'write_data - internal error'
      end select

      if (present(offset)) then
        ! Write Header
        write(text1,'(i0)') ncomps
        write(text2,'(i0)') offset
        write(fid) '    <DataArray Name="', trim(label), '"', &
        &  ' type=', trim(data_text), ' NumberOfComponents="', trim(text1), &
        &  '" format="appended" offset="', trim(text2), '" />', LF
        offset = offset + HEADERTYPE_SIZE + n*nbytes*ncomps

      else
        ! Write Data
        write(fid) int(n*nbytes*ncomps, kind=HEADERTYPE_KIND)
        select case(real_or_int)
        case(META_IS_REAL)
          select case(data_kind) ! real(...) kind must be known at compile time
          case(SP)
            write(fid) real(rdata, kind=SP)
          case(DP)
            write(fid) real(rdata, kind=DP)
          case default
            error stop 'write_data - invalid branch 1'
          end select
        case(META_IS_INT)
          select case(data_kind) ! int(...) kind must be known at compile time
          case(I1B)
            write(fid) int(idata, kind=I1B)
          case(I4B)
            write(fid) int(idata, kind=I4B)
          case(I8B)
            write(fid) int(idata, kind=I8B)
          case default
            error stop 'write_data - invalid branch 2'
          end select
        end select
      end if
    end subroutine write_data


    subroutine extract_and_write_data(choice_export, graph, fid, vtudata, offset)
      integer(I1B), intent(in) :: choice_export
      class(graph_t), intent(in) :: graph
      integer, intent(in) :: fid
      type(vtuio_data_t), intent(in), optional :: vtudata
      integer, intent(inout), optional :: offset
!
! Copy data to rdata/idata arrays and call write_data
! Must be called twice - see write_data for details
!
      integer, allocatable :: idata(:,:)
      real(dp), allocatable :: rdata(:,:)
      integer :: nitems, i, j, choice_import

      block
        integer :: npoints, ncells
        call get_counts(graph, npoints, ncells)
        select case(choice_export)
        case(META_IS_POINT)
          nitems = npoints
          select type(graph)
          class is (mesh_t)
            ! no data written to point data section
            return
          class default
            ! data marked as point data exported from vertices to point data section
            choice_import = META_IS_POINT
          end select
        case(META_IS_CELL)
          nitems = ncells
          select type(graph)
          class is (mesh_t)
            ! data marked as point data imported from vertices to cell data section
            choice_import = META_IS_POINT
          class default
            ! data marked as cell data imported from edges to cell data section
            choice_import = META_IS_CELL
          end select
        case default
          error stop 'extract and write data - invalid data_export'
        end select
      end block
      if (present(vtudata)) then
        if (allocated(vtudata%meta)) then
          do i=1,size(vtudata%meta)
            associate(m=>vtudata%meta(i))
              if (m%iclass==choice_import+META_IS_REAL) then
                ! real data exported within current context
                call reallocate(rdata, [m%ncomp, nitems])

                if (.not. present(offset)) then
                  select case(choice_import)
                  case(META_IS_POINT)
                    do j=1,nitems
                      rdata(:,j) = graph%vertices(j)%rpar(m%start:m%start+m%ncomp-1)
                    end do
                  case(META_IS_CELL)
                    do j=1,nitems
                      rdata(:,j) = graph%edges(j)%rpar(m%start:m%start+m%ncomp-1)
                    end do
                  end select
                end if
                call write_data(fid, m%nbytes, trim(m%label), rdata=rdata, offset=offset)

              else if (m%iclass==choice_import+META_IS_INT) then
                ! integer data exported within current context
                call reallocate(idata, [m%ncomp, nitems])

                if (.not. present(offset)) then
                  select case(choice_import)
                  case(META_IS_POINT)
                    do j=1,nitems
                      idata(:,j) = graph%vertices(j)%ipar(m%start:m%start+m%ncomp-1)
                    end do
                  case(META_IS_CELL)
                    do j=1,nitems
                      idata(:,j) = graph%edges(j)%ipar(m%start:m%start+m%ncomp-1)
                    end do
                  end select
                end if
                call write_data(fid, m%nbytes, trim(m%label), idata=idata, offset=offset)

              end if
            end associate
          end do
        end if
      end if

    end subroutine extract_and_write_data


    subroutine write_time_value(fid, time)
      integer, intent(in) :: fid
      real(DP), intent(in) :: time

      character(len=20) text1
      write(fid) '<FieldData>'//LF
      write(fid) '  <DataArray type="Float32" Name="TimeValue"&
      & NumberOfTuples="1" format="ascii">'
      write(text1,*) real(time, kind=SP)
      write(fid) '  '//trim(adjustl(text1))//LF
      write(fid) '  </DataArray>'//LF
      write(fid) '</FieldData>'//LF
    end subroutine write_time_value


    ! -----------------------------
    ! Reading the Unstructured Grid
    ! -----------------------------
subroutine vtuio_read2(file, graph, mask, time, vtudata)
  character(len=*), intent(in) :: file
  class(graph_t), intent(inout) :: graph
  integer, intent(in) :: mask(4)
    !! mask(1) - index pointing to "radius" in "vertex%rpar" array
    !! mask(2) - index pointing to the first "position" component in "vertex%rpar" array
    !! mask(3) - index pointing to "type" in "vertex%ipar" array
    !! mask(4) - index pointing to "type" in "edge%ipar" array
  real(DP), intent(out), optional :: time
  type(vtuio_data_t), optional :: vtudata
  type(vtuio_data_t) :: vtudata0

  print *, 'WARNING - using obsolete vtuio_read interface'
  if (present(vtudata)) then
    call vtudata%add_item('radius', mask(MASK_RADIUS), int(META_IS_REAL+META_IS_POINT), 1, 4)
    call vtudata%add_item('type', mask(MASK_POINT_TYPE), int(META_IS_INT+META_IS_POINT), 1, 1)
    call vtudata%add_item('con t', mask(MASK_CELL_TYPE), int(META_IS_INT+META_IS_CELL), 1, 1)
    call vtuio_read1(file, graph, mask(MASK_POSITION), time, vtudata)
  else
    call vtudata0%add_item('radius', mask(MASK_RADIUS), int(META_IS_REAL+META_IS_POINT), 1, 4)
    call vtudata0%add_item('type', mask(MASK_POINT_TYPE), int(META_IS_INT+META_IS_POINT), 1, 1)
    call vtudata0%add_item('con t', mask(MASK_CELL_TYPE), int(META_IS_INT+META_IS_CELL), 1, 1)
    call vtuio_read1(file, graph, mask(MASK_POSITION), time, vtudata0)
  end if
end subroutine

    subroutine vtuio_read1(file, graph, position_id, time, vtudata)
      character(len=*), intent(in) :: file
      type(graph_t), intent(inout) :: graph
      integer, intent(in), optional :: position_id
      type(vtuio_data_t), intent(in), optional :: vtudata
      real(DP), intent(out), optional :: time
!
! Read from VTU file
!
      integer :: fid, npoints, ncells, ios, i, n_meta
      integer :: offset_points, offset_cones, offset_offsets, offset_types
      integer, allocatable :: offset_meta(:)
      integer :: time_pos(2), data_pos(2)
      integer(HEADERTYPE_KIND) :: nblock   ! change KIND if error (!)
      integer(CONNECTIONS_KIND) :: vids(2)
      integer, allocatable :: idata(:)
      real(DP) :: time0
      real(POSITIONS_KIND) :: xloc(3)
      real(DP), allocatable :: rdata(:)
      type(object_t), target :: root
      type(object_t), pointer :: grid, piece, opoints, ocells, point_data, &
          cell_data, data_root
      character(len=1) :: ch
      type(handle_t), allocatable :: points(:)
      type(handle_t) :: cone

      time0 = 0.0

      n_meta = 0
      if (present(vtudata)) then
        if (allocated(vtudata%meta)) n_meta = size(vtudata%meta)
      end if
      allocate(offset_meta(n_meta), source=DATAARRAY_NOT_FOUND)

      ! PART ONE
      ! read and analyze the vtk-tree
      call vtuio_tree_read(file//SUFFIX, root)

      grid => root%findtag('UnstructuredGrid')
      if (.not. associated(grid)) error stop &
          'vtuio_read - tag "UnstructuredGrid" not found'

      piece => grid%findtag('Piece')
      if (.not. associated(piece)) error stop &
          'vtuio_read - tag "Piece" not found'

      ! get "npoints" and "ncells"
      npoints = parse_value(piece, 'NumberOfPoints', 'npoints')
      ncells = parse_value(piece, 'NumberOfCells', 'ncells')

      ! get offsets for points and connections
      opoints => piece%findtag('Points')
      if (.not. associated(opoints)) error stop &
          'vtuio_read - tag "Points" not found'
      offset_points = parse_value(opoints, 'offset', 'offset_points')

      ocells => piece%findtag('Cells')
      if (.not. associated(ocells)) error stop &
          'vtuio_read - tag "Cells" not found'
      offset_cones = parse_value(ocells, &
          'Name', 'connectivity', 'offset', 'offset_cones')
      offset_offsets = parse_value(ocells, &
          'Name', 'offsets', 'offset', 'offset_offsets')
      offset_types = parse_value(ocells, &
          'Name', 'types', 'offset', 'offset_types')

      ! is timevalue present?
      block
        character(len=:), allocatable :: text
        logical :: was_found
        time_pos = 0
        text = grid%findval('Name','TimeValue','format', was_found)
        if (was_found) then
          time_pos = grid%findraw()
        end if
#ifdef DEBUG
        print '("time_pos ",i0,1x,i0)', time_pos
#endif
      end block

      ! get offsets for PointData / CellData in vtudata
      point_data => grid%findtag('PointData')
      cell_data => grid%findtag('CellData')

      block
character(len=:), allocatable :: type
integer :: ncomp

        do i=1, size(offset_meta)
          associate(m=>vtudata%meta(i))
            select case (mod(m%iclass,2))
            case(META_IS_CELL)
              if(associated(cell_data)) then
                call parse_dataarray(cell_data, trim(m%label), type, ncomp, &
                    offset_meta(i))
print *, trim(m%label), ncomp, type, offset_meta(i)
              end if
            case(META_IS_POINT)
              if (associated(point_data)) then
                call parse_dataarray(point_data, trim(m%label), type, ncomp, &
                    offset_meta(i))
print *, trim(m%label), ncomp, type, offset_meta(i)
              end if
            case default
              error stop 'vtuio_read - invalid branch'
            end select
            if (offset_meta(i)==DATAARRAY_NOT_FOUND) &
                print '("WARNING vtuio_read - ",a," not found in file")', trim(m%label)
          end associate
        end do
      end block

      ! where binary data start?
      data_root => root%findtag('AppendedData')
      if (.not. associated(grid)) error stop &
      & 'vtuio_read - tag "AppendedData" not found'
      data_pos = data_root%findraw()
#ifdef DEBUG
      print '("data_pos ",i0,1x,i0)', data_pos
#endif


      ! PART TWO
      ! Verify consistency of offsets / npoints / ncells with actual datablock
      open(newunit=fid, file=file//SUFFIX, status='old', access='stream')

      ! Reposition "data_pos" to start at "_" marker
      do
        read(fid, pos=data_pos(1)) ch
        if (ch==' ' .or. ch==LF) then
          data_pos(1) = data_pos(1) + 1
          cycle
        end if
        if (ch /= '_') error stop 'vtuio_read - could not find appended data'
        exit
      end do

      ! First value in each data block is a header indicating the number of
      ! bytes for the datablock. Header can be 4 or 8 bytes integer.
      !
      ! 8 bytes headers read into 4 byte "nblock" will pass validation, but
      ! imported data will be shifted by 1 byte (!!!)
      ! TODO better validation???
      read(fid, pos=data_pos(1)+1+offset_points) nblock
      associate(item=>int(nblock)/(3*npoints), check=>mod(int(nblock),3*npoints))
#ifdef DEBUG
        print '("-points = ",i0,1x,i0,1x,i0)', nblock, item, check
#endif
        if (check/=0) error stop &
            'vtuio_read - validation fails, header size 32/64 mismatch?'
        if (item/=POSITIONS_SIZE) error stop &
            'vtuio_read - validation fails, position kind mismatch'
      end associate

      read(fid, pos=data_pos(1)+1+offset_cones) nblock
      associate(item=>int(nblock)/(2*ncells), check=>mod(int(nblock),2*ncells))
#ifdef DEBUG
        print '("-cones = ",i0,1x,i0,1x,i0)', nblock, item, check
#endif
        if (check/=0) error stop &
            'vtuio_read - validation fails, header size 32/64 mismatch?'
        if (item/=CONNECTIONS_SIZE) error stop &
            'vtuio_read - validation fails, connections kind mismatch'
      end associate

      read(fid, pos=data_pos(1)+1+offset_offsets) nblock
#ifdef DEBUG
      print '("-offsets = ",i0,1x,i0,1x,i0)', nblock, int(nblock)/(ncells), mod(int(nblock),ncells)
#endif

      read(fid, pos=data_pos(1)+1+offset_types) nblock
#ifdef DEBUG
      print '("-types = ",i0,1x,i0,1x,i0)', nblock, int(nblock)/(ncells), mod(int(nblock),ncells)
#endif


      ! PART THREE
      ! Initialize and read points from the file
      call graph%initialize(npoints, ncells, is_directed_graph=.false.)
      allocate(points(npoints))
      read(fid, pos=data_pos(1)+1+offset_points) nblock ! skip header
      block
        type(vertex_t) :: vdummy
        integer :: ipar(size(vdummy%ipar))
        real(DP) :: rpar(size(vdummy%rpar))
        ipar = -77 ! arbitrary values (for debugging)
        rpar = 0.11e-20_dp ! arbitrary values
        do i=1, npoints
          read(fid) xloc
          rpar(position_id:position_id+2) = real(xloc, kind=DP)
          points(i) = graph%add_vertex(ipar, rpar)
        end do
      end block

!     if (offset_r >= 0) then
!       call read_data(fid, data_pos(1)+1+offset_r, npoints, rdata=rdata)
!       graph%vertices(1:npoints)%rpar(mask(MASK_RADIUS)) = rdata
!     else
!       print '("VTU file Warning: radius component is not present")'
!     end if

!     if (offset_at >= 0) then
!       call read_data(fid, data_pos(1)+1+offset_at, npoints, idata=idata)
!       graph%vertices(1:npoints)%ipar(mask(MASK_POINT_TYPE)) = idata
!     else
!       print '("VTU file Warning: point type component is not present")'
!     end if


      ! PART FOUR
      ! Read connections from the file
      read(fid, pos=data_pos(1)+1+offset_cones) nblock ! skip header
      block
        type(edge_t) :: edummy
        integer :: ipar(size(edummy%ipar))
        real(DP) :: rpar(size(edummy%rpar))
        ipar = -42         ! arbitrary values
        rpar = 0.11e-20_dp ! arbitrary values
        do i=1, ncells
          read(fid) vids
          vids = vids + 1 ! our indices start at 1 (not 0)
          if (any(vids<1) .or. any(vids>npoints)) &
            & error stop 'vtuio_read - connection points at non-existent atom'
          cone = graph%add_edge(points(vids(1)), points(vids(2)), ipar, rpar)
        end do
      end block

!     if (offset_ct >= 0) then
!       call read_data(fid, data_pos(1)+1+offset_ct, ncells, idata=idata)
!       graph%edges(1:ncells)%ipar(mask(MASK_CELL_TYPE)) = idata
!     else
!       print '("VTU file Warning: cell type component is not present")'
!     end if

      ! PART FIVE
      ! Read time component
      block
        character(len=:), allocatable :: val
        if (all(time_pos>0)) then
          if (allocated(val)) deallocate(val)
          allocate(character(len=time_pos(2)-time_pos(1)) :: val)
          read(fid, pos=time_pos(1)) val
          read(val,*,iostat=ios) time0
          if (ios/=0) error stop 'vtuio_read - error reading time value'
        end if
      end block

      ! THE END
      print '("VTU file read: points=",i0," cells=",i0)', npoints, ncells
      close(fid)

      if (present(time)) then
        time = time0
        print '("VTU file read: time=",f8.2)', time
      end if
    end subroutine vtuio_read1


    subroutine parse_dataarray(obj, name, type, ncomp, offset)
      type(object_t), intent(in), target :: obj
      character(len=*), intent(in) :: name
      character(len=:), allocatable, intent(out) :: type
      integer, intent(out) :: ncomp, offset
!
! If DataArray could not be found - return offset with -1
! If DataArray was found, but
!   - type attribute not found
!     -> warn and return type as a zero length character
!   - NumberOfComponents attribute not found
!     -> warn and return ncomp with 1
!   - offset attribute not found - data attachment method not supported
!     -> error
!
      character(len=MAX_BUFFER_LEN) :: text
      type(object_t), pointer :: dataarray
      logical :: was_found
      integer :: ios

      if (trim(name)=='') then
        dataarray => obj%findtag('DataArray')
      else
        text = obj%findval('Name', trim(name), 'type', was_found, dataarray)
!TODO may report non-existing DataArray if type not present. This may be misleading
! TODO need more convenient function, not interessed in type/text here
      end if

      if (.not. associated(dataarray)) then
        allocate(character(len=0)::type)
        ncomp = 1
        offset = DATAARRAY_NOT_FOUND
        print '("WARNING - DataArray ",a," not found")', name
        return
      end if

      text = dataarray%findval('type', was_found)
      if (was_found) then
        allocate(character(len=len_trim(text)) :: type)
        type = trim(text)
      else
        allocate(character(len=0) :: type)
        print '("WARNING - type attribute not found in DataArray ",a)', name
      end if

      text = dataarray%findval('NumberOfComponents', was_found)
      if (was_found) then
        read(text,*,iostat=ios) ncomp
        if (ios/=0) error stop &
            'parse_dataarray - invalid number in NumberOfComponents attribute'
      else
        print '("WARNING - NumberOfComponents not present in ",a,&
            & ". Assuming 1")', name
        ncomp = 1
      end if

      text = dataarray%findval('offset', was_found)
      ios = 0
      if (was_found) read(text,*,iostat=ios) offset
      if (ios/=0 .or. .not. was_found) &
        error stop &
          'parse_dataarray - a compulsory offset could not be find/read'
    end subroutine parse_dataarray


    subroutine read_data(fid, pos_start, nvals, rdata, idata)
      integer, intent(in) :: fid, pos_start, nvals
      real(DP), intent(out), allocatable, optional :: rdata(:)
      integer, intent(out), allocatable, optional :: idata(:)
!
! Read binary values into rdata/idata array. "pos_start" is the header
! position. The header value read from file and the number of values "nvals"
! are used to determine whether values are written in 1/4/8 bytes format.
!
      integer, parameter :: IS_REAL=1, IS_INT=2
      integer :: ritype ! IS_REAL or IS_INT
      integer(HEADERTYPE_KIND) :: nblock
      integer :: nbytes
      integer(I1B), allocatable :: idata8(:)
      integer(I4B), allocatable :: idata32(:)
      integer(I8B), allocatable :: idata64(:)
      real(SP), allocatable :: rdata32(:)
      real(DP), allocatable :: rdata64(:)

      if (present(rdata) .and. .not. present(idata)) then
        ritype = IS_REAL
        allocate(rdata(nvals))
      else if (.not. present(rdata) .and. present(idata)) then
        ritype = IS_INT
        allocate(idata(nvals))
      else
        error stop 'read_data - rdata/idata must be present, but not both'
      end if

      read(fid, pos=pos_start) nblock
      nbytes = int(nblock)/nvals
#ifdef DEBUG
      print '("read_data - ",i0," blocks, ",i0," bytes. Zero check ",i0)', &
          nblock, nbytes, mod(int(nblock),nvals)
#endif
      if (mod(int(nblock),nvals)/=0) error stop &
      & 'read_data - validation fails (input data inconsistent)'

      select case(ritype)
      case(IS_REAL)
        select case(nbytes)
        case(4)
          allocate(rdata32(nvals))
          read(fid) rdata32
          rdata = rdata32
        case(8)
          allocate(rdata64(nvals))
          read(fid) rdata64
          rdata = real(rdata64)
        case default
          error stop 'read_data - unsupported real data kind'
        end select
      case(IS_INT)
        select case(nbytes)
        case(1)
          allocate(idata8(nvals))
          read(fid) idata8
          idata = int(idata8)
        case(4)
          allocate(idata32(nvals))
          read(fid) idata32
          idata = idata32
        case(8)
          allocate(idata64(nvals))
          read(fid) idata64
          idata = int(idata64)
        case default
          error stop 'read_data - unsupported integer data kind'
        end select
      end select
    end subroutine read_data


!
! Parse numerical value from the object
!
!TODO move to vtuio_tree?
    function parse_value1(obj, name, msg) result(ivalue)
      type(object_t), intent(in) :: obj
      character(len=*), intent(in) :: name, msg
      integer :: ivalue
      integer :: ios
      logical :: was_found
      character(len=:), allocatable :: text

      ios = -1
      text = obj%findval(name, was_found)
      if (was_found) read(text,*,iostat=ios) ivalue
      if (ios/=0) error stop 'parse_value - could not determine '//msg
#ifdef DEBUG
      print '(a," = ",i0)', msg, ivalue
#endif
    end function parse_value1


    function parse_value2(obj, name1, value1, name2, msg) result(ivalue)
      type(object_t), intent(in) :: obj
      character(len=*), intent(in) :: name1, value1, name2, msg
      integer :: ivalue

      integer :: ios
      logical :: was_found
      character(len=:), allocatable :: text

      ios = -1
      text = obj%findval(name1, value1, name2, was_found)
      if (was_found) read(text,*,iostat=ios) ivalue
      if (ios/=0) error stop 'parse_value - could not determine '//msg
#ifdef DEBUG
      print '(a," = ",i0)', msg, ivalue
#endif
    end function parse_value2


    subroutine reallocate_real(arr, required_shape)
      real(DP), intent(inout), allocatable :: arr(:,:)
      integer, intent(in) :: required_shape(2)

      if (allocated(arr)) then
        if (any(shape(arr) /= required_shape)) deallocate(arr)
      end if
      if (.not. allocated(arr)) &
          allocate(arr(required_shape(1), required_shape(2)))
    end subroutine reallocate_real


    subroutine reallocate_int(arr, required_shape)
      integer, intent(inout), allocatable :: arr(:,:)
      integer, intent(in) :: required_shape(2)

      if (allocated(arr)) then
        if (any(shape(arr) /= required_shape)) deallocate(arr)
      end if
      if (.not. allocated(arr)) &
          allocate(arr(required_shape(1), required_shape(2)))
    end subroutine reallocate_int


    ! ==============================
    ! Organize real and integer data
    ! ==============================

    subroutine meta_add_item(this, label, start, iclass, ncomp, nbytes)
      class(vtuio_data_t), intent(inout) :: this
      character(len=*), intent(in) :: label
      integer, intent(in) :: iclass
        !! use VTUIO_META_POINT/CELL + VTUIO_META_R/I to select the correct
        !! value
      integer, intent(in) :: ncomp, nbytes
      integer, intent(in) :: start

      if (.not. allocated(this%meta)) allocate(this%meta(0))

      ! validate iclass/ncomp
      if (iclass < 0 .or. iclass > 3) &
          error stop 'meta_add_item - invalid iclass'
      if (ncomp /= 1 .and. ncomp /=3) &
          print '("WARNING: expected scalar or 3d-vector")'
      block
        character(len=10) :: data_text
        integer :: data_kind, ubound
        ! error stops if nbytes/iclass combination is invalid
        call data_descriptors(2*(iclass/2), nbytes, data_kind, data_text)
        ! validate start
        select case(iclass)
        case(META_IS_REAL + META_IS_POINT)
          ubound = VSIZE_RPAR
        case(META_IS_REAL + META_IS_CELL)
          ubound = ESIZE_RPAR
        case(META_IS_INT + META_IS_POINT)
          ubound = VSIZE_IPAR
        case(META_IS_INT + META_IS_CELL)
          ubound = ESIZE_IPAR
        case default
          error stop 'meta_add_item - internal error'
        end select
        if (start < 1 .or. start+ncomp-1 > ubound) error stop &
            'meta_add_item - start / iclass combination out of bounds'
      end block

      this%meta = [this%meta, &
          vtuio_meta_t(start,ncomp,nbytes,label,int(iclass,I1B))]
    end subroutine meta_add_item


    subroutine meta_free(this)
      class(vtuio_data_t), intent(inout) :: this
      if (allocated(this%meta)) deallocate(this%meta)
    end subroutine meta_free

  end module vtuio_mod
