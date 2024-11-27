program main
  use omp_lib
  integer(kind=4) :: coreid, group1, group2
  call sleep(5)
  !$omp target 
  !$omp parallel num_threads(4) private(coreid, group1, group2)
  call getcoreid(coreid)
  group1 = omp_get_thread_num()
  print *, "level:", 1, "core:", coreid, "tnum:", omp_get_thread_num()
  !$omp parallel num_threads(4) private(coreid) shared(group1) private(group2)
  call getcoreid(coreid)
  group2 = omp_get_thread_num()
  print *, "level:", 2, "core:", coreid, "tnum:", group1, omp_get_thread_num()
  !$omp parallel num_threads(4) private(coreid) shared(group1, group2)
  call getcoreid(coreid)
  print *, "level:", 3, "core:", coreid, "tnum:", group1, group2, omp_get_thread_num()
  !$omp end parallel
  !$omp end parallel
  !$omp end parallel
  !$omp end target
end program main