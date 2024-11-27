module m
  type :: dt
    integer :: a
  end type dt
  contains
subroutine test(c, a, b, defg)
  integer, optional :: a(:)
  integer, optional, pointer :: b(:)
  integer :: c(:)
  type(dt), optional :: defg
  !$omp target
  print *, present(a), present(b), c, defg%a
  if (present(a)) then
    print *, a
  endif
  !$omp end target
  ! endif
end subroutine test
end module m
program  main
  use m
  implicit none
  integer :: c(3)
  call test(c)
end program  main