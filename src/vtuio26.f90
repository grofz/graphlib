! CONTENT (only public procedures are listed)
! module vtuio_mod
!   subroutine vtuio_write(file, graph, mask, time, vtudata)
!   subroutine vtuio_read(file, graph, mask, time)
!
! THIS IS A MODIFIED VERSION (July 2026)

  module vtuio_mod
    use graph_mod, only : graph_t, handle_t
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

    integer(I1B), parameter, public :: VTUIO_META_I=0, VTUIO_META_R=2,&
        VTUIO_META_POINT=0, VTUIO_META_CELL=1

    type, public :: vtuio_data_t
      !! 1. Call `add_item` to add additional data fields
!TODO not implemented yet
!     !! 2. Call `reallocate` to allocate arrays
!     !! 3. Copy data to `??dat` arrays
      !! 4. Call vtuio_write
      !! 5. Call `finalize`
      integer, allocatable  :: pidat(:,:)
        !! shape = [no of components of data, no of points]
      integer, allocatable  :: cidat(:,:)
        !! shape = [no of components of data, no of cells]
      real(DP), allocatable :: prdat(:,:)
      real(DP), allocatable :: crdat(:,:)
      type(vtuio_meta_t), allocatable :: meta(:)
      integer :: totcomp(0:3) = 0
    contains
      procedure :: add_item => meta_add_item
      procedure :: reallocate => meta_reallocate
      procedure :: free => meta_free
    end type vtuio_data_t

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
! TYPE OF CELLS (all cells are VTK_LINE "3")
    integer, parameter :: CELLTYPE_KIND=I1B, CELLTYPE_SIZE=1
    character(len=*), parameter :: CELLTYPE_TEXT='"UInt8"'
    integer, parameter :: VTK_LINE = 3
! -----------------------------------------------------------------------------

    integer, parameter :: MAX_BUFFER_LEN=400
    character(len=1), parameter :: LF=char(10) ! end-of-line

    character(len=*), parameter :: SUFFIX = '.vtu'

    integer, parameter :: RADIUS_DATA_SIZE = 4

    integer, parameter :: DEBUG=0

    ! legend for items in "mask" array argument
    integer, parameter :: MASK_RADIUS=1, MASK_POSITION=2, MASK_POINT_TYPE=3, MASK_CELL_TYPE=4

  contains

    ! -----------------------------
    ! Writing the Unstructured Grid
    ! -----------------------------

    subroutine vtuio_write(file, graph, mask, time, vtudata)
      !* Write Unstructured Grid - vertices and edges
      character(len=*), intent(in) :: file
        !! file name without .vtu suffix
      type(graph_t), intent(in) :: graph
        !! contains vertices a.k.a. points and edges a.k.a. cells
      integer, intent(in) :: mask(4)
        !! mask(1) - index pointing to "radius" in "vertex%rpar" array
        !! mask(2) - index pointing to the first "position" component in "vertex%rpar" array
        !! mask(3) - index pointing to "type" in "vertex%ipar" array
        !! mask(4) - index pointing to "type" in "edge%ipar" array
      real(DP), intent(in), optional :: time
      type(vtuio_data_t), optional :: vtudata
        !! data structure for additional data

      character(len=MAX_BUFFER_LEN/2) :: ch1, ch2
      integer :: npoints, ncells, i, fid, offset
      real(DP), allocatable :: rdata(:,:)
      integer, allocatable :: idata(:,:)

      npoints = graph%nvertices
      ncells = graph%nedges

      offset = 0
      open(newunit=fid, file=trim(file)//SUFFIX, status='replace', &
      &    access='stream')

      write(fid) '<VTKFile type="UnstructuredGrid" version="1.0"&
      & byte_order="LittleEndian" ', HEADERTYPE_TEXT, '>', LF
      write(fid) '<UnstructuredGrid>', LF

      if (present(time)) call write_time_value(fid, time)

      write(ch1,'(i0)') npoints
      write(ch2,'(i0)') ncells
      write(fid) '<Piece NumberOfPoints="', trim(ch1), &
      &          '" NumberOfCells="', trim(ch2), '">', LF

      ! Point and Cell data headers
      write(fid) '  <PointData>', LF

      allocate(rdata(1,npoints))
      call write_data(fid, RADIUS_DATA_SIZE, 'radius', rdata=rdata, offset=offset)

      allocate(idata(1,npoints))
      call write_data(fid, 1, 'type', idata=idata, offset=offset)

      if (present(vtudata)) then
        if (allocated(vtudata%meta)) then
          do i=1,size(vtudata%meta)
            associate(m=>vtudata%meta(i))
              if (m%iclass==VTUIO_META_POINT+VTUIO_META_R) then
                call write_data(fid, m%nbytes, trim(m%label), rdata=&
                    vtudata%prdat(m%start:m%start+m%ncomp-1,:), offset=offset)
              else if (m%iclass==VTUIO_META_POINT+VTUIO_META_I) then
                call write_data(fid, m%nbytes, trim(m%label), idata=&
                    vtudata%pidat(m%start:m%start+m%ncomp-1,:), offset=offset)
              end if
            end associate
          end do
        end if
      end if

      write(fid) '  </PointData>', LF

      write(fid) '  <CellData>', LF

      if (allocated(idata)) deallocate(idata)
      allocate(idata(1,ncells))

      call write_data(fid, 1, 'con t', idata=idata, offset=offset)

      if (present(vtudata)) then
        if (allocated(vtudata%meta)) then
          do i=1,size(vtudata%meta)
            associate(m=>vtudata%meta(i))
              if (m%iclass==VTUIO_META_CELL+VTUIO_META_R) then
                call write_data(fid, m%nbytes, trim(m%label), rdata=&
                    vtudata%crdat(m%start:m%start+m%ncomp-1,:), offset=offset)
              else if (m%iclass==VTUIO_META_CELL+VTUIO_META_I) then
                call write_data(fid, m%nbytes, trim(m%label), idata=&
                    vtudata%cidat(m%start:m%start+m%ncomp-1,:), offset=offset)
              end if
            end associate
          end do
        end if
      end if

      write(fid) '  </CellData>', LF

      ! Points - header
      write(fid) '  <Points>', LF
      write(ch1,'(i0)') offset
      write(fid) '    <DataArray type=', POSITIONS_TEXT, &
      & ' NumberOfComponents="3" format="appended" offset="', trim(ch1), &
      & '" />', LF
      offset = offset + HEADERTYPE_SIZE + 3*npoints*POSITIONS_SIZE
      write(fid) '  </Points>', LF

      ! Cells (connections) - header
      write(fid) '  <Cells>', LF
      write(ch1,'(i0)') offset
      write(fid) '    <DataArray type=', CONNECTIONS_TEXT, &
      & ' Name="connectivity" format="appended" offset="', trim(ch1), &
      & '" />', LF
      offset = offset + HEADERTYPE_SIZE + 2*ncells*CONNECTIONS_SIZE
      write(ch1,'(i0)') offset
      write(fid) '    <DataArray type=', CONNECTIONS_TEXT, &
      & ' Name="offsets" format="appended" offset="', trim(ch1), '" />', LF
      offset = offset + HEADERTYPE_SIZE + ncells*CONNECTIONS_SIZE
      write(ch1,'(i0)') offset
      write(fid) '    <DataArray type=', CELLTYPE_TEXT, &
      & ' Name="types" format="appended" offset="', trim(ch1), '" />', LF
      offset = offset + HEADERTYPE_SIZE + ncells*CELLTYPE_SIZE
      write(fid) '  </Cells>', LF

      ! Footer
      write(fid) '</Piece>', LF
      write(fid) '</UnstructuredGrid>', LF
      write(fid) '<AppendedData encoding="raw">', LF
      write(fid) '_' ! data block starts with underscore

      ! Binary data
      if (allocated(rdata)) deallocate(rdata)
      allocate(rdata(1,npoints))
      rdata(1,:) = graph%vertices(1:npoints)%rpar(mask(MASK_RADIUS))
      call write_data(fid, RADIUS_DATA_SIZE, 'r', rdata=rdata)

      if (allocated(idata)) deallocate(idata)
      allocate(idata(1,npoints))
      idata(1,:) = graph%vertices(1:npoints)%ipar(mask(MASK_POINT_TYPE))
      call write_data(fid, 1, 'tp', idata=idata)

      if (present(vtudata)) then
        if (allocated(vtudata%meta)) then
          do i=1,size(vtudata%meta)
            associate(m=>vtudata%meta(i))
              if (m%iclass==VTUIO_META_POINT+VTUIO_META_R) then
                call write_data(fid, m%nbytes, trim(m%label), rdata=&
                    vtudata%prdat(m%start:m%start+m%ncomp-1,:))
              else if (m%iclass==VTUIO_META_POINT+VTUIO_META_I) then
                call write_data(fid, m%nbytes, trim(m%label), idata=&
                    vtudata%pidat(m%start:m%start+m%ncomp-1,:))
              end if
            end associate
          end do
        end if
      end if

!vector data example
!     call write_data(fid, 4, 'v', &
!         rdata=transpose(reshape([atoms(:)%v(1),atoms(:)%v(2),atoms(:)%v(3)],[npoints,3])))

      if (allocated(idata)) deallocate(idata)
      allocate(idata(1,ncells))
      idata(1,:) = graph%edges(1:ncells)%ipar(mask(MASK_CELL_TYPE))
      call write_data(fid, 1, 'con t', idata=idata)

      if (present(vtudata)) then
        if (allocated(vtudata%meta)) then
          do i=1,size(vtudata%meta)
            associate(m=>vtudata%meta(i))
              if (m%iclass==VTUIO_META_CELL+VTUIO_META_R) then
                call write_data(fid, m%nbytes, trim(m%label), rdata=&
                    vtudata%crdat(m%start:m%start+m%ncomp-1,:))
              else if (m%iclass==VTUIO_META_CELL+VTUIO_META_I) then
                call write_data(fid, m%nbytes, trim(m%label), idata=&
                    vtudata%cidat(m%start:m%start+m%ncomp-1,:))
              end if
            end associate
          end do
        end if
      end if

      call write_points(fid, graph, mask(MASK_POSITION))
      call write_connectivity(fid, graph)
      write(fid) int(ncells*CONNECTIONS_SIZE, kind=HEADERTYPE_KIND)
      write(fid) (int(i*2, kind=CONNECTIONS_KIND), i=1,ncells)
      write(fid) int(ncells*CELLTYPE_SIZE, kind=HEADERTYPE_KIND)
      write(fid) (int(VTK_LINE, kind=CELLTYPE_KIND), i=1,ncells)

      write(fid) ' ', LF
      write(fid) '</AppendedData>', LF
      write(fid) '</VTKFile>'
      close(fid)

      print '("VTU write: ",i0," points and ",i0," cells written&
          & to """,a,"""")', npoints, ncells, trim(file)//SUFFIX
    end subroutine vtuio_write


    subroutine write_points(fid, graph, positions_id)
      integer, intent(in) :: fid, positions_id
      type(graph_t), intent(in) :: graph

      integer :: i, npoints

      npoints = graph%nvertices

      write(fid) int(3*npoints*POSITIONS_SIZE, kind=HEADERTYPE_KIND)

      do i=1, npoints
        write(fid) real(graph%vertices(i)%rpar(positions_id:positions_id+2), kind=POSITIONS_KIND)
      end do
    end subroutine write_points


    subroutine write_connectivity(fid, graph)
      integer, intent(in) :: fid
      type(graph_t), intent(in) :: graph

      integer :: i, ncells, vids(2)

      ncells = graph%nedges

      write(fid) int(2*ncells*CONNECTIONS_SIZE, kind=HEADERTYPE_KIND)

      do i=1, ncells
        vids = graph%edges(i)%vertex_indices(graph)
        if (any(vids<1) .or. any(vids>graph%nvertices)) &
        &    error stop 'write_connectivity - a connection to unknown atom'
        ! for Paraview, index starts at zero (not at one)
        write(fid) &
             int(vids(1)-1, kind=CONNECTIONS_KIND), &
             int(vids(2)-1, kind=CONNECTIONS_KIND)
      end do
    end subroutine write_connectivity


    subroutine write_data(fid, nbytes, label, rdata, idata, offset)
      integer, intent(in) :: fid
      integer, intent(in) :: nbytes ! 1, 4, 8
      character(len=*), intent(in) :: label
      real(DP), intent(in), optional :: rdata(:,:)
      integer, intent(in), optional :: idata(:,:)
      integer, intent(inout), optional :: offset
!
! This subroutine must be called twice with same arguments and in a correct
! order. In first run (with offset present), the header is written. In second
! run, called within "AppendedData" section, data is written.
!
! Array "rdata" or "idata" must be present (not both)
! - dimension 1 is component (scalars or vectors)
! - dimension 2 is point / connection
! If "offset" is present then write the header, otherwise write data
!
      integer, parameter :: IS_REAL=1, IS_INT=2
      integer :: ritype ! IS_REAL or IS_INT
      integer :: data_kind, data_size, ncomps, n
      character(len=10) :: data_text
      character(len=MAX_BUFFER_LEN/2) :: ch1, ch2

      if (present(rdata) .and. .not. present(idata)) then
        ritype = IS_REAL
      else if (.not. present(rdata) .and. present(idata)) then
        ritype = IS_INT
      else
        error stop 'write_data - rdata/idata must be present, but not both'
      end if

      select case(ritype)
      case(IS_REAL)
        ncomps = size(rdata, dim=1)
        n = size(rdata, dim=2)
        select case(nbytes)
        case(4)
          data_kind = SP
          data_size = 4
          data_text = '"Float32"'
        case(8)
          data_kind = DP
          data_size = 8
          data_text = '"Float64"'
        case default
          error stop 'write_data - 32 or 64 bytes for real data'
        end select
      case(IS_INT)
        ncomps = size(idata, dim=1)
        n = size(idata, dim=2)
        select case(nbytes)
        case(1)
          data_kind = I1B
          data_size = 1
          data_text = '"UInt8"' ! unsigned integer!
        case(4)
          data_kind = I4B
          data_size = 4
          data_text = '"Int32"'
        case(8)
          ! Note: as "idata" is only 4B, this is not necessary (TODO)
          data_kind = I8B
          data_size = 8
          data_text = '"Int64"'
        case default
          error stop 'write_data - 8 32 64 bytes for int data'
        end select
      end select

      if (present(offset)) then
        ! Write Header
        write(ch1,'(i0)') ncomps
        write(ch2,'(i0)') offset
        write(fid) '    <DataArray Name="', trim(label), '"', &
        &  ' type=', trim(data_text), ' NumberOfComponents="', trim(ch1), &
        &  '" format="appended" offset="', trim(ch2), '" />', LF
        offset = offset + HEADERTYPE_SIZE + n*data_size*ncomps

      else
        ! Write Data
        write(fid) int(n*data_size*ncomps, kind=HEADERTYPE_KIND)
        select case(ritype)
        case(IS_REAL)
          select case(data_kind) ! real(...) kind must be known at compile time
          case(SP)
            write(fid) real(rdata, kind=SP)
          case(DP)
            write(fid) real(rdata, kind=DP)
          case default
            error stop 'write_data - invalid branch 1'
          end select
        case(IS_INT)
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


    subroutine write_time_value(fid, time)
      integer, intent(in) :: fid
      real(DP), intent(in) :: time

      character(len=20) ch1
      write(fid) '<FieldData>'//LF
      write(fid) '  <DataArray type="Float32" Name="TimeValue"&
      & NumberOfTuples="1" format="ascii">'
      write(ch1,*) real(time, kind=SP)
      write(fid) '  '//trim(adjustl(ch1))//LF
      write(fid) '  </DataArray>'//LF
      write(fid) '</FieldData>'//LF
    end subroutine write_time_value



    ! -----------------------------
    ! Reading the Unstructured Grid
    ! -----------------------------

    subroutine vtuio_read(file, graph, mask, time)
      character(len=*), intent(in) :: file
      type(graph_t), intent(inout) :: graph
      integer, intent(in) :: mask(4)
        !! mask(1) - index pointing to "radius" in "vertex%rpar" array
        !! mask(2) - index pointing to the first "position" component in "vertex%rpar" array
        !! mask(3) - index pointing to "type" in "vertex%ipar" array
        !! mask(4) - index pointing to "type" in "edge%ipar" array
      real(DP), intent(out), optional :: time
!
! Read from VTU file
!
      integer :: fid, npoints, ncells, ios, i
      integer :: offset_points, offset_cones, offset_offsets, offset_types
      integer :: offset_r, offset_at, offset_ct
      integer :: time_pos(2), data_pos(2)
      integer(HEADERTYPE_KIND) :: nblock   ! change KIND if error (!)
      integer(CONNECTIONS_KIND) :: vids(2)
      integer, allocatable :: idata(:)
      real(DP) :: time0
      real(POSITIONS_KIND) :: xloc(3)
      real(DP), allocatable :: rdata(:)
      type(object_t), target :: root
      type(object_t), pointer :: obj1, obj2, obj3
      character(len=1) :: ch
      character(len=:), allocatable :: val
      logical :: was_found
      type(handle_t), allocatable :: points(:)
      type(handle_t) :: cone

      time0 = 0.0

      ! PART ONE
      ! read and analyze the vtk-tree
      call vtuio_tree_read(file//SUFFIX, root)

      obj1 => root%findtag('UnstructuredGrid')
      if (.not. associated(obj1)) error stop &
      & 'vtuio_read - tag "UnstructuredGrid" not found'

      obj2 => obj1%findtag('Piece')
      if (.not. associated(obj2)) error stop &
      & 'vtuio_read - tag "Piece" not found'

      ! get "npoints" and "ncells"
      ios = -1
      val = obj2%findval('NumberOfPoints', was_found)
      if (was_found) read(val,*,iostat=ios) npoints
      if (ios/=0) error stop 'vtuio_read - could not determine npoints'
      if (DEBUG>0) print *, 'npoints=',npoints
      ios = -1
      val = obj2%findval('NumberOfCells', was_found)
      if (was_found) read(val,*,iostat=ios) ncells
      if (ios/=0) error stop 'vtuio_read - could not determine ncells'
      if (DEBUG>0) print *, 'ncells=',ncells

      ! get offsets for points and connections
      obj3 => obj2%findtag('Points')
      if (.not. associated(obj3)) error stop &
      & 'vtuio_read - tag "Points" not found'
      ios = -1
      val = obj3%findval('offset', was_found)
      if (was_found) read(val,*,iostat=ios) offset_points
      if (ios/=0) error stop 'vtuio_read - could not determine offset_points'
      if (DEBUG>0) print '("offset points ",i0)', offset_points
      obj3 => obj2%findtag('Cells')
      if (.not. associated(obj3)) error stop &
      & 'vtuio_read - tag "Cells" not found'
      ios = -1
      val = obj3%findval('Name','connectivity','offset', was_found)
      if (was_found) read(val,*,iostat=ios) offset_cones
      if (ios/=0) error stop 'vtuio_read - could not determine offset_cones'
      if (DEBUG>0) print '("offset cones ",i0)', offset_cones
      ios = -1
      val = obj3%findval('Name','offsets','offset', was_found)
      if (was_found) read(val,*,iostat=ios) offset_offsets
      if (ios/=0) error stop 'vtuio_read - could not determine offset_offsets'
      if (DEBUG>0) print '("offset offsets ",i0)', offset_offsets
      ios = -1
      val = obj3%findval('Name','types','offset', was_found)
      if (was_found) read(val,*,iostat=ios) offset_types
      if (ios/=0) error stop 'vtuio_read - could not determine offset_types'
      if (DEBUG>0) print '("offset types ",i0)', offset_types

      ! is timevalue present?
      time_pos = 0
      val = obj1%findval('Name','TimeValue','format', was_found)
      if (was_found) then
        time_pos = obj1%findraw()
      end if
      if (DEBUG>0) print '("time pos ",i0,1x,i0)', time_pos

      ! get offsets for point data
      offset_r = -1
      offset_at = -1
      obj2 => obj1%findtag('PointData')
      if (associated(obj2)) then
        ios = -1
        val = obj2%findval('Name','radius','offset', was_found)
        if (was_found) read(val,*,iostat=ios) offset_r
        ios = -1
        val = obj2%findval('Name','type','offset', was_found)
        if (was_found) read(val,*,iostat=ios) offset_at
      end if

      ! get offsets for cell data
      offset_ct = -1
      obj2 => obj1%findtag('CellData')
      if (associated(obj2)) then
        ios = -1
        val = obj2%findval('Name','con t','offset', was_found)
        if (was_found) read(val,*,iostat=ios) offset_ct
      end if
      if (DEBUG>0) then
        print '("offset radii ",i0)', offset_r
        print '("offset atom type ",i0)', offset_at
        print '("offset cell type ",i0)', offset_ct
      end if

      ! where binary data start?
      obj1 => root%findtag('AppendedData')
      if (.not. associated(obj1)) error stop &
      & 'vtuio_read - tag "AppendedData" not found'
      data_pos = obj1%findraw()
      if (DEBUG>0) print '("data pos ",i0,1x,i0)', data_pos


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
      associate(item=>nblock/(3*npoints), check=>mod(nblock,3*npoints))
        if (DEBUG>0) print '("-points = ",i0,1x,i0,1x,i0)', nblock, item, check
        if (check/=0) error stop &
        & 'vtuio_read - validation fails, header size 32/64 mismatch?'
        if (item/=POSITIONS_SIZE) error stop &
        & 'vtuio_read - validation fails, position kind mismatch'
      end associate

      read(fid, pos=data_pos(1)+1+offset_cones) nblock
      associate(item=>nblock/(2*ncells), check=>mod(nblock,2*ncells))
        if (DEBUG>0) print '("-cones = ",i0,1x,i0,1x,i0)', nblock, item, check
        if (check/=0) error stop &
        & 'vtuio_read - validation fails, header size 32/64 mismatch?'
        if (item/=CONNECTIONS_SIZE) error stop &
        & 'vtuio_read - validation fails, connections kind mismatch'
      end associate

      read(fid, pos=data_pos(1)+1+offset_offsets) nblock
      if (DEBUG>0) print '("-offsets = ",i0,1x,i0,1x,i0)', nblock, nblock/(ncells), mod(nblock,ncells)

      read(fid, pos=data_pos(1)+1+offset_types) nblock
      if (DEBUG>0) print '("-types = ",i0,1x,i0,1x,i0)', nblock, nblock/(ncells), mod(nblock,ncells)


      ! PART THREE
      ! Initialize and read points from the file
      call graph%initialize(npoints, ncells, is_directed_graph=.false.)
      allocate(points(npoints))
      read(fid, pos=data_pos(1)+1+offset_points) nblock ! skip header
      block
        integer, allocatable :: ipar(:)
        real(DP), allocatable :: rpar(:)
        allocate(ipar(graph%niv), source=-77)         ! arbitrary values
        allocate(rpar(graph%nrv), source=0.11e-20_dp) ! arbitrary values
        do i=1, npoints
          read(fid) xloc
          rpar(mask(MASK_POSITION):mask(MASK_POSITION)+2) = real(xloc, kind=DP)
          points(i) = graph%add_vertex(ipar, rpar)
        end do
      end block

      if (offset_r >= 0) then
        call read_data(fid, data_pos(1)+1+offset_r, npoints, rdata=rdata)
        graph%vertices(1:npoints)%rpar(mask(MASK_RADIUS)) = rdata
      else
        print '("VTU file Warning: radius component is not present")'
      end if

      if (offset_at >= 0) then
        call read_data(fid, data_pos(1)+1+offset_at, npoints, idata=idata)
        graph%vertices(1:npoints)%ipar(mask(MASK_POINT_TYPE)) = idata
      else
        print '("VTU file Warning: point type component is not present")'
      end if


      ! PART FOUR
      ! Read connections from the file
      read(fid, pos=data_pos(1)+1+offset_cones) nblock ! skip header
      block
        integer, allocatable :: ipar(:)
        real(DP), allocatable :: rpar(:)
        allocate(ipar(graph%nie), source=-42)         ! arbitrary values
        allocate(rpar(graph%nre), source=0.11e-20_dp) ! arbitrary values
        do i=1, ncells
          read(fid) vids
          vids = vids + 1 ! our indices start at 1 (not 0)
          if (any(vids<1) .or. any(vids>npoints)) &
            & error stop 'vtuio_read - connection points at non-existent atom'
          cone = graph%add_edge(points(vids(1)), points(vids(2)), ipar, rpar)
        end do
      end block

      if (offset_ct >= 0) then
        call read_data(fid, data_pos(1)+1+offset_ct, ncells, idata=idata)
        graph%edges(1:ncells)%ipar(mask(MASK_CELL_TYPE)) = idata
      else
        print '("VTU file Warning: cell type component is not present")'
      end if

      ! PART FIVE
      ! Read time component
      if (all(time_pos>0)) then
        if (allocated(val)) deallocate(val)
        allocate(character(len=time_pos(2)-time_pos(1)) :: val)
        read(fid, pos=time_pos(1)) val
        read(val,*,iostat=ios) time0
        if (ios/=0) error stop 'vtuio_read - error reading time value'
      end if

      ! THE END
      print '("VTU file read: points=",i0," cells=",i0)', npoints, ncells
      close(fid)

      if (present(time)) then
        time = time0
        print '("VTU file read: time=",f8.2)', time
      end if
    end subroutine vtuio_read


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
      nbytes = nblock/nvals
      if (DEBUG>0) print '("read_data - ",i0," blocks, ",i0," bytes")', &
          nblock, nbytes, mod(nblock,nvals)
      if (mod(nblock,nvals)/=0) error stop &
      & 'read_data - validation fails (input data incosistent)'

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
          rdata = rdata64
        case default
          error stop 'read_data - unsupported real data kind'
        end select
      case(IS_INT)
        select case(nbytes)
        case(1)
          allocate(idata8(nvals))
          read(fid) idata8
          idata = idata8
        case(4)
          allocate(idata32(nvals))
          read(fid) idata32
          idata = idata32
        case(8)
          allocate(idata64(nvals))
          read(fid) idata64
          idata = idata64
        case default
          error stop 'read_data - unsupported integer data kind'
        end select
      end select
    end subroutine read_data



    ! ==============================
    ! Organize real and integer data
    ! ==============================

    function meta_add_item(this, label, iclass, ncomp, nbytes) result(start)
      class(vtuio_data_t), intent(inout) :: this
      character(len=*), intent(in) :: label
      integer(I1B), intent(in) :: iclass
        !! use VTUIO_META_POINT/CELL + VTUIO_META_R/I to select the correct
        !! value
      integer, intent(in) :: ncomp, nbytes
      integer :: start

      if (.not. allocated(this%meta)) allocate(this%meta(0))

      if (iclass < 0 .or. iclass > ubound(this%totcomp,1)) &
          error stop 'meta_add_item - invalid iclass'
      if (ncomp /= 1 .and. ncomp /=3) &
          print '("WARNING: expected scalar or 3d-vector")'

      start = this%totcomp(iclass) + 1
      this%totcomp(iclass) = this%totcomp(iclass) + ncomp
      this%meta = [this%meta, &
          vtuio_meta_t(start,ncomp,nbytes,label,int(iclass,I1B))]
    end function meta_add_item


    subroutine meta_reallocate(this, npoints, ncells)
      class(vtuio_data_t), intent(inout) :: this
      integer, intent(in) :: npoints, ncells
      allocate(this%pidat(this%totcomp(0), npoints))
      allocate(this%cidat(this%totcomp(1), ncells))
      allocate(this%prdat(this%totcomp(2), npoints))
      allocate(this%crdat(this%totcomp(3), ncells))
    end subroutine meta_reallocate


    subroutine meta_free(this)
      class(vtuio_data_t), intent(inout) :: this
      if (allocated(this%pidat)) deallocate(this%pidat)
      if (allocated(this%cidat)) deallocate(this%cidat)
      if (allocated(this%prdat)) deallocate(this%prdat)
      if (allocated(this%crdat)) deallocate(this%crdat)
      if (allocated(this%meta)) deallocate(this%meta)
      this%totcomp = 0
    end subroutine meta_free

  end module vtuio_mod
