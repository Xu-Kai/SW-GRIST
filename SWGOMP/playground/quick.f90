subroutine test(a, b, n)
  integer :: n, c
  real(8) :: a(n), b(n)
  !$omp target parallel do map(tofrom: c)
  do i=1,n
    a(i) = a(i) + b(i)
  end do
  !$omp end target parallel do
end subroutine test