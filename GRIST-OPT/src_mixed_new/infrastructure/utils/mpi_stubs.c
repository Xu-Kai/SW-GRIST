#ifdef __sw_slave__
#pragma omp delcare target
void mpi_abort_(){
  while (1);
}
#endif
#include <stdio.h>
void empty(){
  #pragma omp target
  {
    puts("!!!");
    mpi_abort_();
  }
}
