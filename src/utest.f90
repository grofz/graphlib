!
! Module for unit testing. The "utest_t" structure keeps the list of
! assertions that passed/failed. A summary is generated at the end.
!
  module utest_mod
    use iso_fortran_env, only : dp=>real64, sp=>real32
    implicit none
    private

    type, public :: utest_t
      private
      integer :: n = 0
      integer :: npass = 0
      type(assert_line_t), pointer :: first => null()
      type(assert_line_t), pointer :: last => null()
    contains
      generic :: assert => assert_equals_integer, assert_equals_logical, &
          assert_equals_dp, assert_equals_sp, &
          assert_equals_integer_arr, assert_equals_integer_4arr
      generic :: within_tolerance => &
          within_tolerance_scalar, within_tolerance_arr
      procedure :: summarize
      procedure :: all_passed
      final :: utest_finalize
      procedure, private :: assert_equals_integer, assert_equals_logical, &
          assert_equals_dp, assert_equals_sp, &
          assert_equals_integer_arr, assert_equals_integer_4arr
      procedure, private :: within_tolerance_scalar, within_tolerance_arr
    end type utest_t

    interface utest_t
      module procedure utest_construct
    end interface utest_t

    type assert_line_t
      logical :: ispass
      character(len=:), allocatable :: msg1
      character(len=:), allocatable :: msg2
      type(assert_line_t), pointer :: next => null()
    end type assert_line_t

  contains

    type(utest_t) function utest_construct()
!
! A dummy constructor (not actually needed).
!
      utest_construct % n = 0
      utest_construct % npass = 0
      utest_construct % first => null()
      utest_construct % last => null()
    end function utest_construct


    subroutine utest_finalize(utest)
      type(utest_t), intent(inout) :: utest
!
! Deatructor for utest_t variables.
!
      type(assert_line_t), pointer :: dele, next

      next => utest % first
      do
        if (.not. associated(next)) exit
        dele => next
        next => dele % next
        deallocate(dele)
      end do
      utest % n = 0
      utest % npass = 0
      utest % first => null()
      utest % last => null()
    end subroutine utest_finalize


    subroutine assert_equals_integer(this, a, b, msg1)
      class(utest_t), intent(inout) :: this
      integer, intent(in) :: a, b
      character(len=*), intent(in) :: msg1
!
! Add line to "utest" table specifying that "a"=="b"
!
      type(assert_line_t) :: line
      character(len=100) :: stra, strb

      line % msg1 = msg1
      write(stra,*) a
      write(strb,*) b
      stra = adjustl(stra)
      strb = adjustl(strb)
      if (a == b) then
        line % ispass = .true.
        line % msg2 = trim(stra) // ' == ' // trim(strb)
      else
        line % ispass = .false.
        line % msg2 = trim(stra) // ' /= ' // trim(strb)
      end if
      call add_assertline(this, line)
    end subroutine assert_equals_integer


    subroutine assert_equals_integer_4arr(this, a, b, p, q, msg1)
      class(utest_t), intent(inout) :: this
      integer, intent(in) :: a(:), b(:), p(:), q(:)
      character(len=*), intent(in) :: msg1
!
! Add line to "utest" table specifying that "a(:)"=="b(:)"
!
      type(assert_line_t) :: line
      character(len=100) :: stra, strb, strp, strq
      integer :: imatch_ap, imatch_bq, imatch_aq, imatch_bp
      integer, allocatable :: pp(:), qq(:)

      imatch_ap = count_match(a,p)
      imatch_aq = count_match(a,q)
      imatch_bp = count_match(b,p)
      imatch_bq = count_match(b,q)
      if (imatch_ap+imatch_bq > imatch_aq+imatch_bp) then
        ! A-P B-Q
        pp = p
        qq = q
      else
        ! A-Q B-P
        pp = q
        imatch_ap = imatch_aq
        qq = p
        imatch_bq = imatch_bp
      end if

      line % msg1 = msg1
      if (size(a)/=size(pp) .or. size(b)/=size(qq)) then
        write(stra,*) size(a)
        write(strb,*) size(b)
        write(strp,*) size(pp)
        write(strq,*) size(qq)
        line % ispass = .false.
        line % msg2 = &
          'size(a)='//trim(adjustl(stra))//' and size(b)='//trim(adjustl(strb))// &
          ' /= '// &
          'size(p)='//trim(adjustl(strp))//' and size(q)='//trim(adjustl(strq))
      else
        write(stra,*) imatch_ap
        write(strb,*) imatch_bq
        write(strp,*) size(pp)
        write(strq,*) size(qq)
        line % msg2 = 'Matching '// &
          trim(adjustl(stra))//' out of '//trim(adjustl(strp))// &
          ' and '// &
          trim(adjustl(strb))//' out of '//trim(adjustl(strq))// &
          ' items'
        if (imatch_ap==size(a) .and. imatch_bq==size(b)) then
          line % ispass = .true.
        else
          line % ispass = .false.
        end if
      end if
      call add_assertline(this, line)
    end subroutine assert_equals_integer_4arr


    function count_match(a, b) result(imatch)
      integer, intent(in) :: a(:), b(:)
      integer :: imatch

      integer, allocatable :: asorted(:), bsorted(:)

      if (size(a) /= size(b)) then
        imatch = 0
      else
        asorted = a
        bsorted = b
        call sort_int(asorted)
        call sort_int(bsorted)
        imatch = count(asorted==bsorted)
      end if
    end function count_match


    subroutine assert_equals_integer_arr(this, a, b, msg1)
      class(utest_t), intent(inout) :: this
      integer, intent(in) :: a(:), b(:)
      character(len=*), intent(in) :: msg1
!
! Add line to "utest" table specifying that "a(:)"=="b(:)"
!
      type(assert_line_t) :: line
      character(len=100) :: stra, strb
      integer :: imatch

      line % msg1 = msg1
      if (size(a)/=size(b)) then
        write(stra,*) size(a)
        write(strb,*) size(b)
        stra = adjustl(stra)
        strb = adjustl(strb)
        line % ispass = .false.
        line % msg2 = 'size(a)='//trim(stra)//' /= size(b)='//trim(strb)
      else
        imatch = count_match(a, b)
        write(stra,*) imatch
        write(strb,*) size(a)
        stra = adjustl(stra)
        strb = adjustl(strb)
        line % msg2 = 'Matching '//trim(stra)//' out of '//trim(strb)//' items'
        if (imatch == size(a)) then
          line % ispass = .true.
        else
          line % ispass = .false.
        end if
      end if
      call add_assertline(this, line)
    end subroutine assert_equals_integer_arr


    subroutine sort_int(arr)
      integer, intent(inout) :: arr(:)
!
! Insert sort of integer arrays
!
      integer :: i, j, x

      do i=2, size(arr)
        x = arr(i)
        j = i
        do while (j > 1)
          if (arr(j-1) <= x) exit
          arr(j) = arr(j-1)
          j = j - 1
        end do
        arr(j) = x
      end do
    end subroutine sort_int


    subroutine assert_equals_logical(this, a, b, msg1)
      class(utest_t), intent(inout) :: this
      logical, intent(in) :: a, b
      character(len=*), intent(in) :: msg1
!
! Add line to "utest" table specifying that "a" .eqv. "b"
!
      type(assert_line_t) :: line
      character(len=100) :: stra, strb

      line % msg1 = msg1
      write(stra,*) a
      write(strb,*) b
      stra = adjustl(stra)
      strb = adjustl(strb)
      if (a .eqv. b) then
        line % ispass = .true.
        line % msg2 = trim(stra) // ' is ' // trim(strb)
      else
        line % ispass = .false.
        line % msg2 = trim(stra) // ' is not ' // trim(strb)
      end if
      call add_assertline(this, line)
    end subroutine assert_equals_logical


    subroutine assert_equals_sp(this, a, b, msg1)
      class(utest_t), intent(inout) :: this
      real(sp), intent(in) :: a, b
      character(len=*), intent(in) :: msg1
      call assert_equals_dp(this, real(a,dp), real(b,dp), msg1)
    end subroutine assert_equals_sp


    subroutine assert_equals_dp(this, a, b, msg1)
      class(utest_t), intent(inout) :: this
      real(dp), intent(in) :: a, b
      character(len=*), intent(in) :: msg1
!
! Add line to "utest" table specifying that "a" near "b"
!
      type(assert_line_t) :: line
      character(len=100) :: stra, strb

      line % msg1 = msg1
      write(stra,*) a
      write(strb,*) b
      stra = adjustl(stra)
      strb = adjustl(strb)
      if (abs(a-b)<= max(1.0_dp,a,b)*5.0*epsilon(1.0_dp)) then
        line % ispass = .true.
        line % msg2 = trim(stra) // ' is ' // trim(strb)
      else
        line % ispass = .false.
        line % msg2 = trim(stra) // ' is not ' // trim(strb)
      end if
      call add_assertline(this, line)
    end subroutine assert_equals_dp


    subroutine within_tolerance_scalar(this, a, b, tol, msg1)
      class(utest_t), intent(inout) :: this
      real(dp), intent(in) :: a, b, tol
      character(len=*), intent(in) :: msg1
!
! TODO Documentation block
!
      character(len=1), parameter :: LE=char(10)
      type(assert_line_t) :: line
      character(len=100) :: stra, strb

      line % msg1 = msg1
      write(stra,*) a
      write(strb,*) b
      stra = adjustl(stra)
      strb = adjustl(strb)
      if (abs(a-b)/max(abs(a),abs(b),tol) <= tol) then
        line % ispass = .true.
        line % msg2 = LE//trim(stra)//' nearly '//LE//trim(strb)
      else
        line % ispass = .false.
        line % msg2 = LE//trim(stra)//' differ '//LE//trim(strb)
      end if
      call add_assertline(this, line)
    end subroutine within_tolerance_scalar


    subroutine within_tolerance_arr(this, a, b, tol, msg1)
      class(utest_t), intent(inout) :: this
      real(dp), intent(in) :: a(:), b(:), tol
      character(len=*), intent(in) :: msg1
!
! TODO Documentation block
!
      character(len=1), parameter :: LE=char(10)
      type(assert_line_t) :: line
      character(len=100) :: stra, strb

      line % msg1 = msg1
      write(stra,'(SP,*(g0,1x))') a
      write(strb,'(SP,*(g0,1x))') b
      stra = adjustl(stra)
      strb = adjustl(strb)
      if (all(abs(a-b)/max(abs(a),abs(b),tol) <= tol)) then
        line % ispass = .true.
        line % msg2 = LE//'['//trim(stra)//'] nearly'//LE//'['//trim(strb)//']'
      else
        line % ispass = .false.
        line % msg2 = LE//'['//trim(stra)//'] differ'//LE//'['//trim(strb)//']'
      end if
      call add_assertline(this, line)
    end subroutine within_tolerance_arr


    pure logical function all_passed(this)
      class(utest_t), intent(in) :: this
      all_passed =  this % n == this % npass
    end function all_passed


    subroutine summarize(this)
      class(utest_t), intent(in) :: this
!
! Print out the list of assertions stored in "this"
!
      type(assert_line_t), pointer :: line

      if (this % n == 0) then
        print '(a)', 'Object does not contain any assertions'
        print *
        return
      endif

      print '(a)', '--------------'
      print '(a)', ' Test summary '
      print '(a)', '--------------'
      print '(a,i0,a,i0,a)', 'Passed ',this % npass,' out of ',this%n,' assertions:'
      line => this % first
      do
        if (.not. associated(line)) exit
        print '(a,l1,1x,a)',  '  ', line % ispass, &
        &   line % msg1//': '//line % msg2
        line => line % next
      end do
      if (this % n == this % npass) then
        print '(a)', 'All tests passed'
      else
        print '(a)', 'Some tests failed'
      end if
      print *
    end subroutine summarize


    subroutine add_assertline(this, line)
      class(utest_t), intent(inout) :: this
      type(assert_line_t), intent(in) :: line
!
! Add a new line to "utest" table
!
      type(assert_line_t), pointer :: new_line

      allocate(new_line)
      new_line = line
      this % n = this % n + 1
      if (new_line % ispass) this % npass = this % npass + 1
      if (associated(this % first)) then
        this % last % next => new_line
      else
        this % first => new_line
      endif
      this % last => new_line
      new_line % next => null()
    end subroutine add_assertline

  end module utest_mod
