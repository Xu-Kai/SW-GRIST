#ifdef __sw_slave__
#include "swdefs.h"
#include "swgomp.h"
#include "tasktree.h"
#include <slave.h>
__uncached long ws_ptr[64][MAX_NESTING_LEVEL+1];
struct swgomp_work_share {
  /*该数据结构记录的循环如下：*/
  /*#pragma omp for schedule(dynamic, chnk)*/
  /*for (int i = start; i < end; i += incr)*/
  long start, end, incr, chnk;
  long *ws_ptr;
};
__ldm int cur_wslev = -1;
__ldm struct swgomp_work_share ws[MAX_NESTING_LEVEL];
extern int omp_get_nested();
long *resolve_wsptr(){
  if (omp_get_nested()) {
    return &(ws_ptr[swgomp_get_rootpe()][level_info[cur_lev].tlev]);
  } else {
    return &(ws_ptr[swgomp_get_rootpe()][0]);
  }
}
// __ldm int my_chnk;
int GOMP_loop_dynamic_next(long *istart, long *iend) {
  int ichnk = faal(ws[cur_wslev].ws_ptr);
  *istart = ws[cur_wslev].start + ichnk * ws[cur_wslev].chnk * ws[cur_wslev].incr;
  *iend = *istart + ws[cur_wslev].chnk * ws[cur_wslev].incr;
  if (*iend > ws[cur_wslev].end) *iend = ws[cur_wslev].end;
  if (*istart >= ws[cur_wslev].end) return 0;
  // cal_locked_printf("chnk: %d %d %d %d %d, %ld\n", _CGN, _MYID, ichnk, *istart, *iend, ws[cur_wslev].ws_ptr - (long*)ws_ptr);

  return 1;
}
void swgomp_loop_dynamic_init(long start, long end, long incr, long chunk_size){
  cur_wslev ++;
  ws[cur_wslev].start = start;
  ws[cur_wslev].end = end;
  ws[cur_wslev].incr = incr;
  ws[cur_wslev].chnk = chunk_size;
  ws[cur_wslev].ws_ptr = resolve_wsptr();
  if (_MYID == swgomp_get_rootpe()) {
    *ws[cur_wslev].ws_ptr = 0;
  }
  // my_chnk = omp_get_thread_num();
  GOMP_barrier();
  // cal_locked_printf("%ld %ld %ld %ld\n", start, end, incr, chunk_size);
}
int GOMP_loop_dynamic_start(long start, long end, long incr, long chunk_size, long *istart, long *iend){
  swgomp_loop_dynamic_init(start, end, incr, chunk_size);
  return GOMP_loop_dynamic_next(istart, iend);
}
void GOMP_loop_end_nowait(){
  cur_wslev --;
}
void GOMP_loop_end(){
  GOMP_loop_end_nowait();
  GOMP_barrier();
}
struct swgomp_loop_arg {
  void (*fn)(void*);
  void *dat;
  long start, end, incr, chunk_size;
};
void gomp_dynamic_loop_dispatch(struct swgomp_loop_arg *arg){
  swgomp_loop_dynamic_init(arg->start, arg->end, arg->incr, arg->chunk_size);
  arg->fn(arg->dat);
}
void GOMP_parallel_loop_dynamic(void (*fn)(void*), void *dat, int num_threads, long start, long end, long incr, long chunk_size, unsigned flags) {
  struct swgomp_loop_arg arg;
  arg.fn = fn;
  arg.dat = dat;
  arg.incr = incr;
  arg.chunk_size = chunk_size;
  arg.start = start;
  arg.end = end;
  GOMP_parallel((void (*)(void*))gomp_dynamic_loop_dispatch, &arg, num_threads, flags);
}
#endif