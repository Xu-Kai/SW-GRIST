subroutine workshare(n, m, a, b)
  real(kind=8)                             :: a(n, m)
  real(kind=8)                             :: b(n, m)
  !$omp target
  !$omp parallel workshare
  a = a*b
  !$omp end parallel workshare
  !$omp end target
end subroutine
program main
  real(kind=8)                             :: a(100, 100)
  real(kind=8)                             :: b(100, 100)
  call workshare(100, 100, a, b)
  print *, a(100, 100), (b(100, 100))
end program