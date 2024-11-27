#include <omp.h>
#include <stddef.h>
extern void set_ptr_uncached(void **ptr);
extern void set_ptr_cached(void **ptr);
int main(){
  long a = 0;
  printf("%p\n", &a);
  #pragma omp target parallel map(tofrom:a) num_threads(4)
  {
    void *ap = &a;
    set_ptr_uncached(&ap);
    my_printf("cpe: %d %p %p\n", omp_get_thread_num(), &a, ap);
    asm volatile(
    #ifdef __sw_slave__
      "faal $31, 0(%0)\n\t"
    #else
      ""
    #endif
    ::"r"(ap));
  }   
  printf("%d\n", a);
}