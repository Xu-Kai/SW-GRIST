#include <stdio.h>
#include "swdefs.h"
#include "jobserver.h"
#include "swgomp-stack.h"
struct job_arg {
  void (*pc)(void *);
  void *arg;
  enum stack_loc stkloc;
};
#ifdef __sw_slave__
extern volatile __ldm struct job_arg job_arg;// = {NULL, NULL, STK_DO_NOT_MOVE};
extern volatile __ldm long job_rply, job_cntr;// = 0;
#else
extern volatile struct job_arg job_arg;
extern volatile long job_rply, job_cntr;
#endif
#ifdef __sw_slave__
// void jobserver(){
//   asm volatile("rcsr %0, 0": "=r"(_MYID));
//   _MYID &= 63;
//   init_saved_sp();
//   job_rply = 0;
//   while (1){
//     unsafe_wait_value(&job_rply, 1);
//     void (*pc)(void *) = job_arg.pc;
//     void *arg = job_arg.arg;
//     nspawn ++;
//     switch_stack_run(pc, arg, job_arg.stkloc);
//     job_arg.pc = NULL;
//     job_rply = 0;
//   }
// }
#endif
void job_spawn_p2p(int tid, void (*pc)(void*), void *arg, enum stack_loc stkloc){
  volatile struct job_arg *remote_arg = cpe_pointer(tid, cpe_symbol_addr((void*)&job_arg));
  volatile long *remote_rply = cpe_pointer(tid, cpe_symbol_addr((void*)&job_rply));
  #ifdef __sw_slave__
  // if (_MYID == 0) printf("%p %p\n", remote_arg, remote_rply);
  #endif
  // printf("remote_addr: %p\n", remote_arg);
  // fflush(stdout);
  remote_arg->stkloc = stkloc;
  remote_arg->arg = arg;
  remote_arg->pc = pc;
  asm volatile("memb\n\t":::"memory");
  *remote_rply += 1;
  asm volatile("memb\n\t":::"memory");
}
#ifdef __sw_slave__
void job_spawn_bcast(int mode, int mask, void (*pc)(void*), void *arg, enum stack_loc stkloc) {
  job_arg.pc = pc;
  job_arg.stkloc = stkloc;
  job_arg.arg = arg;
  rma_bcast((void*)&job_arg, sizeof(job_arg), &job_rply, mode, mask);
}
#endif
#ifdef __sw_host__
volatile int is_waiting_job = 0;
#endif
void job_wait_p2p(int tid){
  #ifdef __sw_host__
  is_waiting_job = 1;
  #endif
  volatile long *remote_rply = cpe_pointer(tid, cpe_symbol_addr((void*)&job_rply));
  volatile long *remote_cntr = cpe_pointer(tid, cpe_symbol_addr((void*)&job_cntr));
  do {
    // printf("%ld %ld\n", *remote_rply, *remote_cntr);
  } while (*remote_cntr - *remote_rply <= 0);
  // unsafe_wait_value(remote_rply, 0);
  #ifdef __sw_host__
  is_waiting_job = 0;
  #endif
}

#ifdef __sw_host__
#include <signal.h>
void swgomp_show_slave_wait(int signo){
  char hn[4096], msg[4096];
  // if (is_waiting_job) {
  gethostname(hn, 4096);
  sprintf(msg, "%s-CG%d is_waiting_job=%d", hn, arrayid(), is_waiting_job); 
  puts(msg);
  // }
}
void swgomp_set_slave_wait_signal(int n){
  struct sigaction sa;
  sigaction(n, NULL, &sa);
  sa.sa_handler = swgomp_show_slave_wait;
  sa.sa_flags |= SA_NODEFER;
  sigaction(n, &sa, NULL);
}
void swgomp_set_slave_wait_signal_(int *n) {
  swgomp_set_slave_wait_signal(*n);
}
uint64_t athread_idle();
int __real_athread_init();
int CRTS_init();
void __real___real_athread_spawn(void (*pc)(void*), void *arg, int flush);
void hang(){
  while (1);
}
#include <stdlib.h>
int jobserver_initialized = 0;
void jobserver_athread_init(){
  // atexit(hang);
  if (!jobserver_initialized) {
    for (int i = 0; i < 64; i ++) {
      *(int64_t*)(0x800000180080L | ((int64_t)arrayid())<<40L | i<<24L) |= 1L << 54;
    }
    jobserver_initialized = 1;
    shared_stack_init();
    __real_athread_init();
    __real___real_athread_spawn(slave_jobserver, NULL, 0);
    // puts("b4 crts");
    CRTS_init();
    // puts("aft crts");
  }
}
void __wrap_athread_init() __attribute__((alias("jobserver_athread_init")));
#endif
#ifdef __sw_host__
void print_ptr_(void *a){
  printf("mpe: %p\n", a);
}
void print_sp_(){
  void *sp;
  asm volatile("mov $30, %0" : "=r"(sp));
  // if (sp < swgomp_shared_stack)
  // printf("sp of %d is at %p %p\n", _MYID, sp, swgomp_shared_stack);
}
#else
void print_ptr_(void *a){
  cal_locked_printf("%d %p\n", _MYID, a);
}
void print_sp_(){
  void *sp;
  asm volatile("mov $30, %0" : "=r"(sp));
  cal_locked_printf("sp of %d is at %p %p\n", _MYID, sp, swgomp_shared_stack);
}
#endif
#ifdef TEST
// extern int athread_init();
// extern int __real_athread_spawn(void (*pc)(void*), void *arg, int flush);
// extern int athread_join();
extern __ldm int _MYID;
#ifdef __sw_slave__
#include <slave.h>
#include "tasktree.h"
void print_me(int n){
  printf("L3: %d %d %02x %02x\n", _MYID, level_info[cur_lev].ithr, level_info[cur_lev].rmask & 255, level_info[cur_lev].cmask & 255);
}
void spawn_3(){
  printf("L2: %d %d %02x %02x\n", _MYID, level_info[cur_lev].ithr, level_info[cur_lev].rmask & 255, level_info[cur_lev].cmask & 255);
  task_tree_start(4, print_me, NULL);
}
void spawn_2(){
  printf("L1: %d %d %02x %02x\n", _MYID, level_info[cur_lev].ithr, level_info[cur_lev].rmask & 255, level_info[cur_lev].cmask & 255);
  task_tree_start(4, spawn_3, NULL);
}
void spawn_1(){
  task_tree_start(4, spawn_2, NULL);
}
#endif
extern void slave_jobserver();
extern void slave_spawn_1();

#ifdef __sw_host__
extern void athread_init();
extern void CRTS_init();
// #include <athread.h>
int main(){
  puts("initialize...");
  athread_init();
  CRTS_init();
  puts("init done");
  __real_athread_spawn(slave_jobserver, 0, 0);
  puts("mpe write");
  job_spawn_p2p(0, slave_spawn_1, 0);
  puts("mpe wrote");
  job_wait_p2p(0);
}
#endif
#endif