program test1
  use conts_mod, only : stack_t, queue_t
  implicit none(type, external)

  integer, parameter :: CS = 3
  integer :: i

  type(stack_t) :: s
  type(queue_t) :: q

  call s%initialize(chunksize=CS, capacity=3)
  call q%initialize(chunksize=CS, capacity=3)

  call q%enqueue([1,10,100])
  call q%enqueue([2,20,200])
  print *, q%dequeue()
  call q%enqueue([3,30,300])
  call q%enqueue([4,40,400])
  call q%enqueue([5,50,500])
  print *, q%dequeue()
  call q%enqueue([6,60,600])
  call q%enqueue([7,70,700])
  print *, q%dequeue()
  print *, q%dequeue()
  print *, q%dequeue()
  print *, q%dequeue()
  print *, q%dequeue()
  do i=1,10
    call q%enqueue([i,i*10,i*100])
  end do
  do i=1,5
    print *, q%dequeue()
  end do
  do i=11,15
    call q%enqueue([i,i*10,i*100])
  end do
  do i=6,15
    print *, q%dequeue()
  end do
  print *, 'queue size =', q%size()

  end program