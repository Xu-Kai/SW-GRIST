subroutine test(a, b, c)

! io 
   use omp_lib
   real(kind=8), dimension(1000) :: a, b, c
   call print_ptr(a)
   call print_ptr(b)
   call print_ptr(c)
!$omp target 
    call print_ptr(a)
    call print_ptr(b)
    call print_ptr(c)
!$omp parallel private(ie,length_of_triangle,ilev,pv_at_edge,cell_sum) 
!$omp master
    call print_ptr(a)
    call print_ptr(b)
    call print_ptr(c)
!$omp end master
!$omp end parallel
!$omp end target

   return
  end subroutine test
program main
  real (kind=8) :: a(1000, 3), b(1000, 3), c(1000, 3)
  call test(a(1,1), b(1,1), c(1,1))
  call test(a(1,2), b(1,2), c(1,2))
  call test(a(1,3), b(1,3), c(1,3))
end program main