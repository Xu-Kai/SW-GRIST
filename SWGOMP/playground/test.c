extern int cal_locked_printf(const char *s, ...);
extern int getcoreid_(int *);
#include <omp.h>
#ifdef __sw_host__
volatile int dbg = 0;
#else
extern volatile int dbg;
#endif
int main(){
  // sleep(3);
  #pragma omp target parallel for num_threads(4)
  for (int i = 0; i < 128; i ++) {
    int coreid;
    getcoreid_(&coreid);
    long sp;
    __asm__("mov $30, %0\n\t" :"=r"(sp));
    cal_locked_printf("%d %d %d %p\n", coreid, omp_get_thread_num(), i, sp);
  }
}