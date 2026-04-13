program test2
  use conts_mod, only : pqueue_t, handle_t
  implicit none(type,external)

  type(pqueue_t) :: pq
  type(handle_t) :: h1, h2, h3, h4
  integer :: a(3), i

  call pq%initialize(chunksize=3, capacity=3)

  do i=1,2
    h1 = pq%insert([1,10,100],4)
    h2 = pq%insert([2,20,200],2)
    h3 = pq%insert([3,30,300],3)
    h4 = pq%insert([4,40,400],1)
print *, 'Export ', pq%export()
print *, h1, h2, h3, h4
print *, 'Export priorities ', pq%export_priorities()
print *, 'Export handles ', pq%export_handles()
    a = pq%pop()
    print *, a, '444'
    if (i==1) then
      call pq%update_priority(h1, 8)
    else
      call pq%update_priority(h1, 1)
    end if
    a = pq%pop()
    print *, a, '222'
    if (pq%contains(h1)) then
      call pq%update_priority(h1, 1)
    else
      print *, 'not present anymore'
    end if
    print *, 'size ', pq%size()

    call pq%clear()
  end do
end program