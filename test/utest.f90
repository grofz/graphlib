!
! Module for unit testing. The "utest_t" structure keeps the list of
! assertions that passed/failed. A summary is generated at the end.
!
  module utest_mod
    implicit none
    private

    type, public :: utest_t
      private
      integer :: n = 0
      integer :: npass = 0
      type(assert_line_t), pointer :: first => null()
      type(assert_line_t), pointer :: last => null()
    contains
      generic :: assert_eq => assert_equals_integer, assert_equals_logical
      procedure :: summarize
      final :: utest_finalize
      procedure, private :: assert_equals_integer, assert_equals_logical
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
