!
! Parse routines to process input data
!
! inspired also from: https://github.com/jacobwilliams/AoC-2021/blob/master/src/aoc_utilities.f90
!
! Note for myself: never use character(len=:) arrays of characters, likely bugs in older compilers.
!                  array of string_t seems a more robust way
!
  module parse_mod
    use iso_fortran_env, only : I8B => int64
    implicit none
    private
    public string_t
    public read_strings, split, split_nonempty
    public read_numbers, read_pattern
    public unique_sort
    !public read_pattern

    integer, parameter :: DEFAULT_LINE_LENGTH = 80

    type string_t
      character(len=:), allocatable :: str
    contains
      procedure :: to_int => string_to_int
      procedure :: to_int128 => string_to_int128
    end type string_t
    interface string_t
      module procedure string_new
    end interface

  contains

! =======================
! String_t implementation
! =======================
    type(string_t) function string_new(str) result(this)
      character(len=*), intent(in) :: str
      allocate(character(len=len(str)) :: this%str)
      this%str=str
    end function

    elemental integer function string_to_int(this) result(int)
      class(string_t), intent(in) :: this
      if (.not. allocated(this%str)) error stop 'string_to_int - error, string unallocated'
      read(this%str,*) int
    end function string_to_int

    elemental integer(I8B) function string_to_int128(this) result(int)
      class(string_t), intent(in) :: this
      if (.not. allocated(this%str)) error stop 'string_to_int - error, string unallocated'
      read(this%str,*) int
    end function string_to_int128



    subroutine split(str, separator, items)
      character(len=*), intent(in) :: str, separator
      type(string_t), intent(out), allocatable :: items(:)

      integer, parameter :: INIT_SEPS_MAX=50
      integer :: i, j, len_str, len_sep, n_seps, n_seps_max
      integer :: i1, i2
      integer, allocatable :: iseps(:)

      len_sep = len(separator)
      len_str = len(str)

      ! count, how many times the separator appears in the string and save its positions
      n_seps_max = INIT_SEPS_MAX
      allocate(iseps(n_seps_max))
      n_seps = 0
      j = 1
      do
        if (j > len_str) exit
        i = index(str(j:), separator)
        if (i<=0) exit
        if (n_seps == n_seps_max) call expand(iseps, n_seps_max)
        n_seps = n_seps + 1
        iseps(n_seps) = i+j-1
        j = j + i + (len_sep-1)
      end do

      allocate(items(n_seps+1))

      do i = 1, n_seps+1
        if (i==1) then
          i1 = 1
        else
          i1 = iseps(i-1)+len_sep
        end if

        if (i==n_seps+1) then
          i2 = len_str
        else
          i2 = iseps(i)-1
        end if

        if (i2 >= i1) then
          items(i) = string_t(str(i1:i2))
        else
          items(i) = string_t('')
        end if
      end do

    contains
      subroutine expand(arr, nmax)
        integer, allocatable, intent(inout) :: arr(:)
        integer, intent(inout) :: nmax
        integer, allocatable :: wrk(:)
        integer :: nmax_ext
        if (size(arr)/=nmax) error stop 'split/expand - array dimension does not match with input'
        nmax_ext = nmax*2
        allocate(wrk(nmax_ext))
        wrk(1:nmax) = arr
        nmax = nmax_ext
        call move_alloc(wrk, arr)
      end subroutine expand

    end subroutine split



    subroutine split_nonempty(str, separator, items)
      character(len=*), intent(in) :: str, separator
      type(string_t), intent(out), allocatable :: items(:)
!
! Split first. Then remove empty strings
!
      type(string_t), allocatable :: items0(:)
      integer :: i, n, j

      call split(str, separator, items0)
      n = 0
      do i=1, size(items0)
        if (len_trim(items0(i)%str)/=0) n = n + 1
      end do
      allocate(items(n))
      j = 0
      do i=1, size(items0)
        if (len_trim(items0(i)%str)==0) cycle
        j = j + 1
        items(j) = string_t(items0(i)%str)
      end do
if (j /= n) error stop 'split_nonempty - check failed'

    end subroutine split_nonempty



! ======================
! Reading data from file
! ======================

    function read_strings(file) result(strings)
      character(len=*), intent(in) :: file
      type(string_t), allocatable :: strings(:)

      integer :: fid, i, n
      character(len=:), allocatable :: line

      open(newunit=fid, file=file, status='old')
      n = count_lines(fid)
      allocate(strings(n))
      do i=1, n
        call read_line_from_file(fid, line)
        strings(i) = string_t(line)
      end do
      close(fid)
    end function read_strings



    function read_numbers(file) result(a)
      character(len=*), intent(in) :: file
      integer, allocatable :: a(:)
!
! Read from file with a single integer on every line
!
! (just a demo, how use "read_strings" and "str2int" elementals)
!
      type(string_t), allocatable :: lines(:)
      lines = read_strings(file)
      a = lines % To_int()
    end function read_numbers



    function read_pattern(file, trailing_spaces) result(aa)
      character(len=*), intent(in) :: file
      logical, intent(in), optional :: trailing_spaces
      character(len=1), allocatable :: aa(:,:)
!
! Read "character" 2D matrix from the file.
!
      integer :: fid, nrow, ncol, ios, i, nread
      character(len=5000) :: line
      logical :: mode2

      ! Regular or iregular pattersn
      mode2 = .false.
      if (present(trailing_spaces)) mode2 = trailing_spaces

      ! read no of rows, cols, and make sure all rows are of the same
      ! length
      open(newunit=fid, file=file, status='old')
      nrow = 0
      ncol = -1
      do
        if (.not. mode2) then
          ! regular case
          read(fid,'(a)',iostat=ios) line
          if (ios /= 0) exit
          if (ncol==-1) ncol = len_trim(line)
          if (len_trim(line) /= ncol) &
            error stop 'read_pattern - not all lines have same length'
        else
          ! blanks at the end of lines are important
          read(fid,'(a)',advance='no',size=nread,iostat=ios) line
          if (is_iostat_end(ios) .or. is_iostat_eor(ios)) then
            if (nread==0) exit
            if (ncol==-1) ncol = nread
            if (nread /= ncol) &
              error stop 'read_pattern - not all lines have same length'
          else if (ios==0) then
            error stop 'input line too long'
          else
            error stop 'reading error'
          end if
        end if
        nrow = nrow + 1
      end do

      allocate(aa(ncol,nrow))
      rewind(fid)
      do i=1,nrow
        read(fid,'(*(a))') aa(:,i)
      end do
      close(fid)
      aa = transpose(aa)
    end function read_pattern



! =======================
! Unique sort
! =======================

  function unique_sort(vals) result(outs)
    integer, intent(in) :: vals(:)
    integer, allocatable :: outs(:)
!
! Return array of soreted unique values from "vals"
! (source of inspiration: stackoverflow)
!
    integer, allocatable :: outs0(:)
    integer :: min_val, max_val, i

    allocate(outs0(size(vals))) ! maximum possible no of uniques
    i = 0
    min_val = minval(vals)-1
    max_val = maxval(vals)
    do
      if (min_val >= max_val) exit
      i = i + 1
      min_val = minval(vals, mask=vals>min_val)
      outs0(i) = min_val
    end do
    allocate(outs(i), source=outs0(1:i))

  end function unique_sort




 ! =======================
 ! Local helper procedures
 ! =======================



    function count_lines(fid) result(n)
      integer, intent(in) :: fid
      integer :: n
!
! Get number of lines in the file. Empty lines are included in the count.
!
      integer :: ios
      character(len=1) :: line
      rewind(fid)
      n = 0
      do
        read(fid,'(a)',iostat=ios) line
        if (ios /= 0) exit
        n = n + 1
      end do
      rewind(fid)

    end function count_lines



    function count_lines_block(fid) result(narr)
      integer, intent(in) :: fid
      integer, allocatable :: narr(:)
!
! Get number of lines for each block. Blocks are separated by an empty line
!
      integer :: ios, iblock
      character(len=1) :: line
      rewind(fid)
      narr = [0]
      iblock = 1
      do
        read(fid,'(a)',iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line)==0) then
          iblock = iblock + 1
          narr = [narr, 0]
          cycle
        end if
        narr(iblock) = narr(iblock) + 1
      end do
      rewind(fid)

    end function count_lines_block



    subroutine read_line_from_file(fid, line)
      integer, intent(in) :: fid
      character(len=:), allocatable, intent(out) :: line

      integer :: nread, ios
      character(len=DEFAULT_LINE_LENGTH) :: buffer

      nread = 0
      buffer = ''
      line = ''

      do
        read(fid,'(a)',advance='no',size=nread,iostat=ios) buffer
        if (is_iostat_end(ios) .or. is_iostat_eor(ios)) then
          if (nread>0) line = line//buffer(1:nread)
          exit
        else if (ios==0) then
          line = line//buffer
        else
          error stop 'read_line_from_file - read error'
        end if
      end do

    end subroutine read_line_from_file


  end module parse_mod
