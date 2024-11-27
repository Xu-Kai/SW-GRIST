#include <stdlib.h>
#include <stdint.h>
#include "swdefs.h"
#include "jobserver.h"
#include "taskarray.h"
union array_spawn_arg {
  struct {
    int nthr, cpe0;
    void (*fn)(void *);
    void *arg;
    int row_mask, nteam, iteam, ldmstack;
  } s;
  int512 v;
};

#ifdef __sw_host__
extern __ldm volatile union array_spawn_arg array_spawn_arg;
void taskarray_arg_set(int tid, void (*pc)(void*), void *arg, int nteam, int ldmstack){
  volatile union array_spawn_arg *spawn_arg = cpe_pointer(tid, (void*)&array_spawn_arg);
  spawn_arg->s.nteam = nteam;
  spawn_arg->s.fn = pc;
  spawn_arg->s.arg = arg;
  spawn_arg->s.ldmstack = ldmstack;
  asm volatile("memb\n\t":::"memory"); 
}
#endif
#ifdef __sw_slave__
__ldm struct taskarray_info taskarray_info;
__ldm volatile union array_spawn_arg array_spawn_arg;
__ldm volatile long arg_bcast_rply = 0;
__ldm uint64_t team_bitmap;
extern __ldm char _CGN;
// void my_flush_slave_cache(){
//   for (int i = 0; i < 128*1024; i += 32){
//     __asm__ volatile("flushd $31, 0(%0)\n\tmemb\n\t"::"r"(i):"memory");
//   }
// }
void taskarray_spawn_teams(){
  // flush_slave_cache();
  // my_flush_slave_cache();
  // evict_slave_cache_cont(array_spawn_arg.s.arg, array_spawn_arg.s.arg+128);
  // return;
  if (_MYID == 0) {
    // set_ptr_uncached(&array_spawn_arg.s.arg);
    rma_bcast((void*)&array_spawn_arg, sizeof(array_spawn_arg), &arg_bcast_rply, RMA_MODE_COL_BCAST, 0x1);
    job_spawn_bcast(RMA_MODE_COL_BCAST, 0x1, taskarray_spawn_teams, NULL, STK_IN_SHR);
    // cal_locked_printf("spawn arg: %p\n", array_spawn_arg.s.arg);
    if (_CGN == 0) {
      long *arg = array_spawn_arg.s.arg;
      for (int i = 0; i < 16; i ++) {
        // cal_locked_printf("scarg[%d]: %p %p\n", i, arg[i], &arg[i]);
      }
      // set_ptr_uncached(&arg);
      // for (int i = 0; i < 16; i ++) {
      //   cal_locked_printf("sucarg[%d]: %p\n", i, arg[i]);
      // }
    }
  }
  
  int nteam = array_spawn_arg.s.nteam;
  if (nteam > 8) 
    nteam = 8;
  if (nteam < 1)
    nteam = 1;
  int p2 = 31 - __builtin_clz(nteam);//1:0; 2:1; 4:2; 8:3...
  if (nteam > (1 << p2)) p2++;
  int log2_team_nrow = 3 - p2;
  int chkmask = (1 << log2_team_nrow) - 1;
  taskarray_info.iteam = _ROW >> log2_team_nrow;
  taskarray_info.nteam = nteam;
  int iam_team_head = (_ROW & chkmask) == 0;
  //1: ff, 2: f, 4: 3
  taskarray_info.row_mask = (1<<(1 << log2_team_nrow)) - 1;
  // cal_locked_printf("%d %d\n", _MYID, team_id);
  asm volatile("sync %0\n\t"::"r"(0xff) :"memory");
  team_bitmap = 0;
  for (int i = 0; i < nteam; i ++){
    team_bitmap |= 1UL << ((i<<log2_team_nrow)<<3);
  }
  if (iam_team_head && taskarray_info.iteam < nteam) {
    asm volatile("sync %0\n\t"::"r"(0xff) :"memory");
    if (_CGN == 0) {
      // cal_locked_printf("spawn arg: %d %d %d %x %x %p\n", nteam, taskarray_info.iteam, log2_team_nrow, chkmask, taskarray_info.row_mask, array_spawn_arg.s.arg);
      long *arg = array_spawn_arg.s.arg;
      for (int i = 0; i < 32; i ++) {
        // cal_locked_printf("scarg[%d]: %p\n", i, arg[i]);
      }
    }
    array_spawn_arg.s.fn(array_spawn_arg.s.arg);
    int syn_mask = 0;
    for (int i = 0; i < nteam; i ++){
      syn_mask |= 1 << (i << log2_team_nrow);
    }
    flush_slave_cache();
    asm volatile("sync %0\n\t"::"r"(syn_mask) :"memory");
  } else {
    asm volatile("sync %0\n\t"::"r"(0xff) :"memory");
  }
}
void taskarray_quick_spawn(){
  if (_MYID == 0) {
    rma_bcast((void*)&array_spawn_arg, sizeof(array_spawn_arg), &arg_bcast_rply, RMA_MODE_COL_BCAST, 0xff);
    job_spawn_bcast(RMA_MODE_COL_BCAST, 0xff, taskarray_quick_spawn, NULL, STK_IN_SHR);
  }
  asm volatile("sync %0\n\tsynr %0"::"r"(0xff) :"memory");
  int nthr = array_spawn_arg.s.nteam;
  taskarray_info.iteam = 0;
  taskarray_info.nteam = 1;
  taskarray_info.row_mask = 0xff;
  team_bitmap = -1L;
  taskarray_info.nthr = nthr;
  taskarray_info.ithr = _MYID;
  asm volatile("sync %0\n\tsynr %0"::"r"(0xff) :"memory");

  if (_MYID < nthr) {
    array_spawn_arg.s.fn(array_spawn_arg.s.arg);
    flush_slave_cache();
    asm volatile("sync %0\n\tsynr %0"::"r"(0xff) :"memory");
  } else {
    asm volatile("sync %0\n\tsynr %0"::"r"(0xff) :"memory");
  }
}
void task_array_thread_run(){
  void *dat = array_spawn_arg.s.arg;
  void (*fn)(void*) = array_spawn_arg.s.fn;
  taskarray_info.row_mask = array_spawn_arg.s.row_mask;
  flush_slave_cache();
  task_array_sync();

  int ithr = _MYID - array_spawn_arg.s.cpe0;
  taskarray_info.ithr = ithr;
  taskarray_info.nthr = array_spawn_arg.s.nthr;
  taskarray_info.cpe0 = array_spawn_arg.s.cpe0;
  taskarray_info.nteam = array_spawn_arg.s.nteam;
  taskarray_info.iteam = array_spawn_arg.s.iteam;
  int ldmstack = array_spawn_arg.s.ldmstack;
  if (ithr < array_spawn_arg.s.nthr)
    switch_stack_run(fn, dat, ldmstack ? STK_IN_LDM : STK_IN_SHR);
  flush_slave_cache();
  task_array_sync();
  // asm volatile("sync %0\n\t"
  //              "synr %1\n\t"
  //              ::"r"(taskarray_info.row_mask), "r"(0xff) : "memory");
}
extern void flush_slave_cache();
void task_array_spawn_threads(int nthr, void (*fn)(void*), void *dat){
  // my_flushd_();
  if (_CGN == 0 && _MYID == 0){
    // cal_locked_printf("subarg: %p\n", dat);
    long *arg = dat;
    for (int i = 0; i < 16; i ++) {
      // cal_locked_printf("scarg[%d]: %p\n", i, arg[i]);
    }
  }
  array_spawn_arg.s.nthr = nthr;
  array_spawn_arg.s.cpe0 = _MYID;
  array_spawn_arg.s.arg = dat;
  array_spawn_arg.s.fn = fn;
  array_spawn_arg.s.row_mask = taskarray_info.row_mask;
  // cal_locked_printf("spawn row: %x %p\n", array_spawn_arg.s.row_mask, dat);
  array_spawn_arg.s.iteam = taskarray_info.iteam;
  array_spawn_arg.s.nteam = taskarray_info.nteam;
  rma_bcast((void*)&array_spawn_arg, sizeof(array_spawn_arg), &arg_bcast_rply, RMA_MODE_ROW_BCAST, taskarray_info.row_mask);
  job_spawn_bcast(RMA_MODE_ROW_BCAST, taskarray_info.row_mask, task_array_thread_run, NULL, STK_IN_SHR);
  // flush_slave_cache();
  task_array_thread_run();
}
void task_array_sync(){
  // cal_locked_printf("synm: %x\n", taskarray_info.row_mask);
  asm volatile("sync %0\n\t"
               "synr %1\n\t"
               ::"r"(taskarray_info.row_mask), "r"(0xff) : "memory");
}
#endif