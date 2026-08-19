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


  module vtuio_tree_mod
    !* To parse VTK file and extract the information needed to load data from
    !  binary VTU file.
    !
    ! Version 1 (__.03.2023)
    !
    use iso_fortran_env, only: IOSTAT_END
    implicit none (type, external)
    private
    public object_t, vtuio_tree_read

    integer, parameter :: IDEBUG = 0

integer :: counter=0

    integer, parameter :: MAXSTRLEN=20

    integer, parameter :: TP_RAW=1, TP_AV=2, TP_OBJ=3
    !* To parse VTU file, there are three object types
    !
    !  1. raw data as a range of positions in a file
    !
    !  2. attribute="value" object
    !
    !  3. compound (list of objects of 1/2/3 types)

    type object_t
      private
      integer :: tp=0
        !! type of object: raw data, attribute-value, compound
      integer :: pos(2)=-1
        !! for raw data: range of positions in a file
      character(len=MAXSTRLEN) :: nam='NULL', val='NULL'
        !! for attribute-value: name and value of the attribute
      character(len=MAXSTRLEN) :: tag='NULL'
        !! for compound: tag identification
      type(object_ptr), allocatable :: list(:)
        !! for compound: array of child objects
    contains
      final :: object_finalize
      procedure :: findtag => find_tag
      procedure, private :: find_val1, find_val2
      generic :: findval => find_val1, find_val2
      procedure :: findraw => find_raw
      procedure :: print=>object_print
    end type object_t
    interface object_t
      module procedure new_raw, new_av, new_obj
    end interface

    type object_ptr
      private
      type(object_t), pointer :: ptr
    end type object_ptr

  contains

    ! -----------------------------------------
    ! Tree structure constructors and finalizer
    ! -----------------------------------------

    function new_raw(pos) result(new)
      !! Allocate new raw-data object
      integer, intent(in) :: pos(2)
      type(object_t), pointer :: new

      allocate(new)
counter = counter+1
      new%tp = TP_RAW
      new%pos = pos
    end function new_raw


    function new_av(nam, val) result(new)
      !! Allocate new attribute-value object
      character(len=*), intent(in) :: nam, val
      type(object_t), pointer :: new

      allocate(new)
counter = counter+1
      new%tp = TP_AV
      new%nam = nam
      new%val = val
      if (len_trim(nam)>MAXSTRLEN .or. len_trim(val)>MAXSTRLEN) &
      & error stop 'new_av - name or value are too long, increse MAXSTRLEN'
    end function new_av


    function new_obj(tag) result(new)
      !! Allocate new compound object
      character(len=*), intent(in) :: tag
      type(object_t), pointer :: new

      allocate(new)
counter = counter+1
      new%tp = TP_OBJ
      new%tag = tag
      if (len_trim(tag)>MAXSTRLEN) &
      & error stop 'new_obj - tag too long, increase MAXSTRLEN'
      allocate(new%list(0))
    end function new_obj


    subroutine addobj(this, added)
      !! Add object to a compound object
      class(object_t), intent(inout) :: this
      type(object_t), pointer :: added

      type(object_ptr) :: dat

      if (this%tp /= TP_OBJ) error stop 'addobj - only TP_OBJ can be used'
      if (.not. allocated(this%list)) error stop 'addobj - list not allocated'
      dat%ptr => added
      this%list = [this%list, dat]
    end subroutine addobj


    recursive subroutine object_finalize(this)
      type(object_t), intent(inout) :: this

      integer :: i
      type(object_t), pointer :: obj

      if (this%tp/=TP_OBJ) return
      if (.not. allocated(this%list)) return
      do i=1, size(this%list)
        !deallocate(this%list(i)%ptr) ! internal compiler error...
        obj => this%list(i)%ptr       ! ... is a working alternative
        deallocate(obj)
counter = counter - 1
      end do
      ! must deallocate the list to avoid a "second" finalization (!)
      deallocate(this%list)
      if (IDEBUG>0) print *, 'Finalized tag ',this%tag, counter
    end subroutine


    ! ----------------
    ! Parsing VTU file
    ! ----------------

    subroutine vtuio_tree_read(file, root)
      !! Top-level procedure for parsing the file
      character(len=*), intent(in) :: file
      type(object_t), intent(out) :: root

      integer :: fid, pos, maxpos, icyc
      type(object_t), pointer :: obj

      root%tp = TP_OBJ
      root%tag = 'ROOT'
      allocate(root%list(0))
      open(newunit=fid, file=file, status='old', access='stream', &
      &    form='unformatted')

      do icyc=1,1
        call skip_xmlheader(fid)
        call read_obj(fid, obj)
        call addobj(root, obj)
        call skip_spaces(fid, pos)
        inquire(fid, size=maxpos)
        if (pos > maxpos) exit
      end do
      if (icyc>1) &
      &  print *, 'vtuiotree_read WARNING finished, but not at the end of file'

      close(fid)
print '("File size ",i0," bytes. Scan completed.")', pos-1
    end subroutine vtuio_tree_read


    recursive subroutine read_obj(fid, new)
      !! Construct TP_OBJ object (tag and list of children objects)
      integer, intent(in) :: fid
      type(object_t), pointer, intent(out) :: new

      character(len=MAXSTRLEN) :: tag, nam, val, chtag
      character(len=1) :: lastch
      logical :: was_success, was_closed, err
      type(object_t), pointer :: obj
      integer :: pos

      call read_tag(fid, tag, lastch)
      new => object_t(tag)
      if (IDEBUG>0) print '("Reading tag ",a)', '"'//trim(new%tag)//'"...'

      was_closed = .false.
      if (lastch /= '>') then  ! i.e. not <Tag> but <Tag ... >
        ! Read and add attributes
        do
          call read_attribute(fid, nam, val, was_success, was_closed)
          if (.not. was_success) exit
          call addobj(new, object_t(nam, val))
          if (IDEBUG>0) call object_print(new%list(size(new%list))%ptr)
        end do
      end if

      ! Read children objects
      ! there are three possibilities
      ! - raw data
      ! - closing tag, e.g. </tag>
      ! - another object, e.g. <child-object>
      do
        if (was_closed) exit
        inquire(fid, pos=pos)
        call read_tag(fid, chtag, lastch, err)
        if (err) then
          call read_raw(fid, new%tag, obj)
          call addobj(new, obj)
          if (IDEBUG>0) call object_print(new%list(size(new%list))%ptr)
        else if (chtag(1:1)=='/' .and. chtag(2:)==tag) then
          was_closed=.true.
        else
          read(fid, pos=pos)
          call read_obj(fid, obj) ! recursion (!)
          call addobj(new, obj)
        end if
      end do
      if (IDEBUG>0) call object_print(new)
      if (IDEBUG>0) print '("..reading done ",a,1x,i0)','"'//trim(new%tag)//'"', size(new%list)
    end subroutine read_obj


    subroutine read_attribute(fid, nam, val, was_success, was_closed)
      integer, intent(in) :: fid
      character(len=*), intent(out) :: nam, val
      logical, intent(out) :: was_success, was_closed

      character(len=1) :: ch
      integer :: n

      call skip_spaces(fid)
      read(fid) ch

      ! return if no more attributes
      if (ch=='/') then
        read(fid) ch
        if (ch /= '>') error stop 'read_attribute - expected > after /'
        was_success = .false.
        was_closed = .true.
        return
      else if (ch=='>') then
        was_success = .false.
        was_closed = .false.
        return
      end if

      ! read attribute name (until =)
      n = 0
      nam = ''
      do
        if (ch=='=') exit
        n = n+1
        if (n>len(nam)) &
        &  error stop 'read_attribute - name too long, increase MAXSTRLEN'
        nam(n:n) = ch
        read(fid) ch
      end do

      call skip_spaces(fid)
      read(fid) ch

      ! read attribute value (closed in "...")
      if (ch /= '"') error stop 'read_attribute - " expected'
      n = 0
      val = ''
      do
        read(fid) ch
        if (ch=='"') exit
        n = n+1
        if (n>len(val)) &
        &  error stop 'read_attribute - value too long, increase MAXSTRLEN'
        val(n:n) = ch
      end do

      was_success = .true.
    end subroutine read_attribute


    subroutine read_tag(fid, tag, lastch, err)
      integer, intent(in) :: fid
      character(len=*), intent(out) :: tag
      character(len=1), intent(out) :: lastch
      logical, intent(out), optional :: err

      integer :: n, pos
      character(len=1) :: ch

      inquire(fid, pos=pos)
      call skip_spaces(fid)
      read(fid) ch
      if (ch /= '<') then
        if (present(err)) then
          read(fid, pos=pos)
          err = .true.
          return
        else
          print *, 'ERROR: expected "<", read "'//ch//'"'
          error stop 'read_tag - unexpected character read'
        end if
      end if
      tag = ''
      n = 0
      do
        read(fid) ch
        ! markers for the end of tag: space, LF, ">"
        if (ch==' ' .or. ch==char(10) .or. ch==">") exit
        n = n+1
        if (n>len(tag)) &
        & error stop 'read_tag - tag too long, increase MAXSTRLEN'
        tag(n:n) = ch
      end do
      lastch = ch
      if (present(err)) err = .false.
    end subroutine read_tag


    subroutine read_raw(fid, tag, obj)
      !! Create raw object and skip file until closing "tag" has been found.
      integer, intent(in) :: fid
      character(len=*), intent(in) :: tag
      type(object_t), pointer, intent(out) :: obj

      integer :: pos0, pos1
      character(len=len_trim(tag)+3) :: ctag

      inquire(fid, pos=pos0)
      ctag = '</'//trim(tag)//'>'
      call scan_pattern(fid, ctag)
      inquire(fid, pos=pos1)
      obj => object_t([pos0, pos1])
    end subroutine read_raw


    subroutine scan_pattern(fid, pattern)
      !! Scan file until a pattern has been found.
      !! Position file before pattern.
      integer, intent(in) :: fid
      character(len=*), intent(in) :: pattern

      integer :: pos, n, ios
      character(len=1) :: ch

      inquire(unit=fid, pos=pos)
      n = 1
      do
        read(fid, iostat=ios) ch
        if (ios==IOSTAT_END) then
          print *, 'ERROR: pattern "'//pattern//'" not found'
          error stop 'scan_pattern - end of file reached'
        end if
        if (ch==pattern(n:n)) then
          n = n + 1
        else
          inquire(unit=fid, pos=pos)
          n = 1
        end if
        if (n==len(pattern)+1) exit
      end do
      read(fid, pos=pos)
    end subroutine scan_pattern


    subroutine skip_xmlheader(fid)
      !! Skip "<?xml ... ?>" at the start of file (if present).
      integer, intent(in) :: fid

      character(len=*), parameter :: xml='<?xml'
      character(len=len(xml)) :: cha
      character(len=*), parameter :: pattern='?>'
      character(len=len(pattern)) :: chx
      integer :: pos

      call skip_spaces(fid)
      inquire(fid, pos=pos)
      read(fid) cha
      if (cha==xml) then
        call scan_pattern(fid, pattern)
        read(fid) chx
print *, 'Xml header skipped'
      else
print *, 'No Xml header'
        read(fid, pos=pos)
      end if
    end subroutine skip_xmlheader


    subroutine skip_spaces(fid, pos)
      integer, intent(in) :: fid
      integer, intent(out), optional :: pos
!
! Skip spaces and LF until another character (or end-of-file) is reached.
!
      character(len=1) :: ch
      integer :: ios, pos0

      do
        inquire(unit=fid, pos=pos0)
        read(fid, iostat=ios) ch
        if (ios == IOSTAT_END) exit
        if (ch==' ' .or. ch==char(10)) cycle
        read(fid, pos=pos0)
        exit
      end do
      if (present(pos)) pos = pos0
    end subroutine skip_spaces


    ! ------------------------------------
    ! Extracting information from the tree
    ! ------------------------------------

    subroutine object_print(this)
      class(object_t), intent(in) :: this

      select case(this%tp)
      case(TP_OBJ)
        print '("  Tag ",a," with ",i0," objects")', trim(this%tag), size(this%list)
      case(TP_AV)
        print '(2x,a," = ",a)', trim(this%nam), trim(this%val)
      case(TP_RAW)
        print '(2x, "Raw start at ",i0," size ",i0)', this%pos(1), this%pos(2)-this%pos(1)
      case default
        error stop 'invalid tp'
      end select
    end subroutine object_print


    recursive function find_tag(this, tag) result(obj)
      class(object_t), intent(in), target :: this
      character(len=*), intent(in) :: tag
      type(object_t), pointer :: obj
!
! Return pointer to a first object with matching "tag" or null pointer.
!
      integer :: i

      obj => null()
      if (this%tp /= TP_OBJ) return
      if (this%tag == tag) then
        obj => this
        return
      end if
      do i=1, size(this%list)
        associate(chl=>this%list(i)%ptr)
          if (chl%tp /= TP_OBJ) cycle
          if (chl%tag == tag) then
            obj => this%list(i)%ptr
          else
            obj => find_tag(chl, tag)
          end if
          if (associated(obj)) exit
        end associate
      end do
    end function find_tag


    recursive function find_val1(this, nam1, val1, nam2, was_found, &
          obj_found) result(val2)
      class(object_t), intent(in), target :: this
      character(len=*), intent(in) :: nam1, nam2, val1
      logical, intent(out) :: was_found
      type(object_t), pointer, intent(out), optional :: obj_found
      character(len=:), allocatable :: val2
!
! Return value of attribute named "nam2" for a first object that has
! attribute "nam1" with value "val1"
!
      integer :: i
      logical :: val1_match, nam2_found
      character(len=:), allocatable :: val2_a, val2_b

      was_found = .false.
      if (present(obj_found)) obj_found => null()
      val1_match = .false.
      nam2_found = .false.
      if (this%tp /= TP_OBJ) return
      do i=1,size(this%list)
        associate(chl=>this%list(i)%ptr)
          if (chl%tp==TP_AV) then
            if (chl%nam==nam1 .and. chl%val==val1) val1_match = .true.
            if (chl%nam==nam2) then
              val2_a = chl%val
              nam2_found = .true.
            end if
          else if (chl%tp==TP_OBJ) then
            val2_b = find_val1(chl, nam1, val1, nam2, was_found, obj_found)
            if (was_found) val2 = val2_b
          end if
          if (was_found) exit
        end associate
      end do
      if (.not. was_found) then
        if (val1_match .and. nam2_found) then
          val2 = val2_a
          was_found = .true.
          if (present(obj_found)) obj_found => this
        end if
      end if
    end function find_val1


    recursive function find_val2(this, nam2, was_found) result(val2)
      class(object_t), intent(in) :: this
      character(len=*), intent(in) :: nam2
      logical, intent(out) :: was_found
      character(len=:), allocatable :: val2
!
! Return value of attribute named "nam2" for a first object.
!
      integer :: i

      was_found = .false.
      if (this%tp /= TP_OBJ) return
      do i = 1, size(this%list)
        associate(chl=>this%list(i)%ptr)
          if (chl%tp==TP_AV) then
            if (chl%nam==nam2) then
              val2 = chl%val
              was_found = .true.
            end if
          else if (chl%tp==TP_OBJ) then
            val2 = find_val2(chl, nam2, was_found)
          end if
        end associate
        if (was_found) exit
      end do
    end function find_val2


    recursive function find_raw(this) result(pos)
      class(object_t), intent(in) :: this
      integer :: pos(2)
!
! Return range of first raw data block found or "-1" if no raw data
!
      integer :: i

      pos = -1
      if (this%tp == TP_RAW) then
        pos = this%pos
        return
      else if (this%tp /= TP_OBJ) then
        return
      end if

      do i=1, size(this%list)
        associate(chl=>this%list(i)%ptr)
          pos = find_raw(chl)
        end associate
        if (all(pos>0)) exit
      end do
    end function find_raw

  end module vtuio_tree_mod
