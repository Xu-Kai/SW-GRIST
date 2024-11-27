#include <stdio.h>
#include <omp.h>

#ifdef __sw_slave__
#include <slave.h>
extern __ldm int _MYID;
#else
extern int _MYID;
#endif
int main(){
  #pragma omp target parallel num_threads(4) 
  {
    int coreid;
    asm volatile("rcsr %0, 0\n\t" : "=r"(coreid));
    _MYID = coreid;
    int group1 = omp_get_thread_num();
    my_printf("level1: core: %p %3d, tunm: %3d\n", &_MYID, _MYID, group1);
    #pragma omp parallel num_threads(4) private(coreid)
    {
      asm volatile("rcsr %0, 0\n\t" : "=r"(coreid));
      int group2 = omp_get_thread_num();
      my_printf("level2: core: %p %3d, tunm: %3d %3d\n", &_MYID, _MYID, group1, group2);
      #pragma omp parallel num_threads(4) private(coreid)
      {
        asm volatile("rcsr %0, 0\n\t" : "=r"(coreid));
        my_printf("level3: core: %p %3d, tunm: %3d %3d %3d\n", &_MYID, _MYID, group1, group2, omp_get_thread_num());
      }
    }
  }
}