#include <stdint.h>
#include <stdlib.h>
#include "swdefs.h"
#include "swgomp-stack.h"
#ifdef __sw_host__
extern __uncached void *swgomp_shared_stack;
extern void *swgomp_shared_stack_alloc;
extern size_t swgomp_stksz;
#endif
#ifdef __sw_slave__
__uncached void *swgomp_shared_stack = NULL;
void *swgomp_shared_stack_alloc = NULL;
size_t swgomp_stksz;
__ldm void *ssps[STK_DO_NOT_MOVE];
#endif
#ifdef __sw_host__
static uint64_t parse_stksz() {
  uint64_t omp_stksz = 8 * 1024L * 1024;
  const char *env_stksz = getenv("OMP_STACKSIZE");
  if (env_stksz) {
    uint64_t omp_stksz_parse = 0;
    int parse_fail = 0;
    for (const char *c = env_stksz; *c; c++) {
      if (*c >= '0' && *c <= '9')
        omp_stksz_parse = omp_stksz_parse * 10 + *c - '0';
      else {
        if (c[1] == '\0') {
          switch (c[0]) {
          case 'B':
            break;
          case 'K':
            omp_stksz_parse <<= 10;
            break;
          case 'M':
            omp_stksz_parse <<= 20;
            break;
          case 'G':
            omp_stksz_parse <<= 30;
            break;
          default:
            parse_fail = 1;
          }
        } else {
          parse_fail = 1;
        }
      }
    }
    if (!parse_fail) {
      omp_stksz = omp_stksz_parse;
    }
  }
  // printf("SWGOMP: stack size is %d\n", omp_stksz);
  return omp_stksz;
}
extern void slave_init_saved_sp();
void shared_stack_init(){
  // if (swgomp_shared_stack == NULL){
    swgomp_stksz = parse_stksz();
    swgomp_shared_stack_alloc = malloc(swgomp_stksz*64+128);
    swgomp_shared_stack = (void*)(((uint64_t)swgomp_shared_stack_alloc) + 63UL & ~63UL);
    // printf("%p\n", swgomp_shared_stack);
    // if (swgomp_shared_stack_alloc == NULL) {
      // puts("!!!");
    // }
  // }
}
#endif
void *shared_stack_get(int pe) {
  // my_printf("%p %p\n", swgomp_shared_stack, (void*)((uint64_t)swgomp_shared_stack + (pe+1) * swgomp_stksz));
  return (void*)((uint64_t)swgomp_shared_stack + (pe+1) * swgomp_stksz);
}
#ifdef __sw_slave__

static enum stack_loc stack_where(){
  void *sp;
  void *const ldm_top = (void*)0x40000;
  asm volatile("mov $30, %0\n\t":"=r"(sp));
  if (sp <= ldm_top) {
    return STK_IN_LDM;
  } else if (sp > shared_stack_get(_MYID-1) || sp <= shared_stack_get(_MYID)) {
    return STK_IN_SHR;
  } else {
    return STK_IN_PRIV;
  }
}
extern long get_slave_cache_size();
extern long get_slave_ldmshare_size();
void init_saved_sp(){
  if (stack_where() != STK_IN_SHR) {
    // cal_locked_printf("SHR %p\n", shared_stack_get(_MYID));
    ssps[STK_IN_SHR] = shared_stack_get(_MYID);
  }
  if (stack_where() != STK_IN_LDM){
    ssps[STK_IN_LDM] = (void*)0x40000 - get_slave_cache_size() - get_slave_ldmshare_size();
    // cal_locked_printf("LDM %p\n", (void*)0x40000 - get_slave_cache_size() - get_slave_ldmshare_size());
  }
  if (stack_where() != STK_IN_PRIV)
    ssps[STK_IN_PRIV] = NULL;
  // for (int i = 0; i < STK_DO_NOT_MOVE; i ++) cal_locked_printf("%d %d %p\n", i, stack_where(), ssps[i]);
}
void switch_stack_run(void (*pc)(void *), void *arg, enum stack_loc loc) {
  enum stack_loc oldloc = stack_where();
  // cal_locked_printf("%d %p %p %d\n", _MYID, pc, arg, loc);
  if (!pc) {
    while (1);
  }
  if (oldloc == loc || loc == STK_DO_NOT_MOVE) {
    pc(arg);
  } else {
    set_stack_run(pc, arg, ssps[loc], ssps + oldloc);
  }
}
#endif