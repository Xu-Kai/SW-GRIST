#include "swdefs.h"
#include "jobserver.h"
struct athread_spawn_arg{
  void (*pc)(void*);
  void *arg;
  int flush;
};
extern void evict_slave_cache_cont(void *a, void *b);
#ifdef __sw_slave__
__ldm struct athread_spawn_arg athread_spawn_arg;
__ldm long athread_spawn_rply = 0;
void athread_spawn_handler(struct athread_spawn_arg *argp){
  if (_MYID == 0) {
    // evict_slave_cache_cont(argp, argp+1);
    flush_slave_cache();
    athread_spawn_arg.pc = argp->pc;
    athread_spawn_arg.arg = argp->arg;
    athread_spawn_arg.flush = argp->flush;
    rma_bcast(&athread_spawn_arg, sizeof(athread_spawn_arg), &athread_spawn_rply, RMA_MODE_ROW_BCAST, 0xff);
    job_spawn_bcast(RMA_MODE_ROW_BCAST, 0xff, (void (*)(void*))athread_spawn_handler, argp, STK_DO_NOT_MOVE);
  } else {
    unsafe_wait_value(&athread_spawn_rply, 1);
    athread_spawn_rply = 0;
  }
  // cal_locked_printf("%p\n", athread_spawn_arg.pc);
  // __asm__ volatile("sync %0\n\tsynr %0\n\t"::"r"(0xff):"memory");
  athread_spawn_arg.pc(athread_spawn_arg.arg);
  if (athread_spawn_arg.flush) flush_slave_cache();
  __asm__ volatile("sync %0\n\tsynr %0\n\t"::"r"(0xff):"memory");
}
#endif
#ifdef __sw_host__
static struct athread_spawn_arg harg;
extern void slave_athread_spawn_handler(void *);
void jobserver_athread_spawn(void (*pc)(void*), void *arg, int flush){
  harg.pc = pc;
  harg.arg = arg;
  harg.flush = flush;
  job_spawn_p2p(0, slave_athread_spawn_handler, &harg, STK_DO_NOT_MOVE);
}
void jobserver_athread_join(){
  job_wait_p2p(0);
}
void __wrap___real_athread_spawn(void (*pc)(void*), void *arg, int flush) __attribute__((alias("jobserver_athread_spawn")));
void __wrap_athread_join() __attribute__((alias("jobserver_athread_join")));
#endif