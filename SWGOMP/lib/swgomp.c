#include "tasktree.h"
#include "jobserver.h"
#include "taskarray.h"
#include "swgomp-stack.h"
#include "gomp-constants.h"

int __real_athread_spawn(void *func, void *arg, int flush_slave_dcache);
#ifdef __sw_slave__
__ldm int swgomp_nested = 0;
void GOMP_parallel(void (*fn)(void*), void *arg, unsigned nthr, unsigned flags) {
  void *argp = arg;
  if ((intptr_t)argp < 0x3FFFF) {
    argp = cpe_pointer(_MYID, argp);
  }
  if (nthr == 0) {
    nthr = 64;
  }

  if (swgomp_nested) {
    if (cur_lev != -1 && nthr > (2 << level_info[cur_lev].next_lev)) {
      nthr = 2 << level_info[cur_lev].next_lev;
    }
    task_tree_start(nthr, fn, argp);
  } else {
    if (nthr == (unsigned)-1) {
      fn(arg);
    } else {
      task_array_spawn_threads(nthr, fn, argp);
    }
  }
}
void getcoreid_(int *x) {
  *x = _MYID;
}
int swgomp_get_cpe_num(){
  return _MYID;
}
int swgomp_get_rootpe(){
  return swgomp_nested ? (level_info[cur_lev].cpe0) : (taskarray_info.cpe0);
}
int omp_get_thread_num(){
  return swgomp_nested ? (cur_lev >= 0 ? level_info[cur_lev].ithr : 0) : taskarray_info.ithr;
}
int omp_get_thread_num_() __attribute__((alias("omp_get_thread_num")));
int omp_get_num_threads(){
  return swgomp_nested ? (cur_lev >= 0 ? level_info[cur_lev].nthr : 1) : taskarray_info.nthr;
}
int omp_get_num_threads_() __attribute__((alias("omp_get_num_threads")));
void omp_set_nested(int nested){
  swgomp_nested = !!nested;
}
int omp_get_nested(){
  return swgomp_nested;
}
int omp_get_team_num(){
  return taskarray_info.iteam;
}
void GOMP_teams (unsigned int num_teams, unsigned int thread_limit){

}
void GOMP_barrier(){
  if (swgomp_nested) {
    task_tree_sync();
  } else {
    task_array_sync();
  }
}

#endif
#ifdef __sw_host__
extern void slave_taskarray_spawn_teams();
extern void slave_taskarray_quick_spawn();
int swgomp_initialized = 0;
void (*SWGOMP_pre_target_cb)(int dev, void *fn, int b, void *dat, void *sz,
                        void *kinds, int c, size_t d, void **args) = NULL;
void (*SWGOMP_post_target_cb)(int dev, void *fn, int b, void *dat, void *sz,
                        void *kinds, int c, size_t d, void **args) = NULL;
void SWGOMP_target_ext(int dev, void *fn, int b, void *dat, void *sz,
                        void *kinds, int c, size_t d, void **args){
  jobserver_athread_init();
  if (SWGOMP_pre_target_cb) SWGOMP_pre_target_cb(dev, fn, b, dat, sz, kinds, c, d, args);
  int num_teams = 0;
  int thread_limit;
  while (*args){
    intptr_t id = (intptr_t) *args ++, val;
    if (id & GOMP_TARGET_ARG_SUBSEQUENT_PARAM) 
      val = (intptr_t) *args ++;
    else
      val = id >> GOMP_TARGET_ARG_VALUE_SHIFT;
    id &= GOMP_TARGET_ARG_ID_MASK;
    if (id == GOMP_TARGET_ARG_NUM_TEAMS)
      num_teams = val;
    else if (id == GOMP_TARGET_ARG_THREAD_LIMIT)
      thread_limit = val;
  }

  if (num_teams == -1) {
    taskarray_arg_set(0, fn, dat, thread_limit, dev != 1);
    job_spawn_p2p(0, slave_taskarray_quick_spawn, NULL, STK_IN_SHR);
    job_wait_p2p(0);
  } else if (num_teams == 0) {
    job_spawn_p2p(0, fn, dat, STK_IN_SHR);
    job_wait_p2p(0);
  } else {
    long *arg = dat;
    // printf("target arg: %p\n", dat);
    for (int i = 0; i < 16; i ++) {
      // if (arrayid() == 0) printf("arg[%d]: %p %p\n", i, arg[i], &arg[i]);
    }

    taskarray_arg_set(0, fn, dat, num_teams, dev != 1);
    job_spawn_p2p(0, slave_taskarray_spawn_teams, NULL, STK_IN_SHR);
    job_wait_p2p(0);
    // printf("target arg: %p\n", dat);
    for (int i = 0; i < 16; i ++) {
      // if (arrayid() == 0) printf("arg[%d]: %p %p\n", i, arg[i], &arg[i]);
    }
  }
  if (SWGOMP_post_target_cb) SWGOMP_post_target_cb(dev, fn, b, dat, sz, kinds, c, d, args);
}
void __wrap_GOMP_target_ext(int a, void *fn, int b, void *dat, void *sz,
                        void *kinds, int c, size_t d, void **args) __attribute__((alias("SWGOMP_target_ext")));
#endif