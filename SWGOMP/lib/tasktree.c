#ifdef __sw_slave__
#include "swdefs.h"
#include "jobserver.h"
#include "tasktree.h"
extern void omp_set_nested(int);
extern int omp_get_nested();
extern void flush_slave_cache();
union tree_spawn_arg {
  struct {
    uint64_t bitmap;
    int clev, tlev;
    int nthr;
    void (*pc)(void *);
    void *arg;
    int nested;
  } s;
  int512 v;
};
__ldm union tree_spawn_arg tasktree_spawn_arg;
__ldm struct level_info level_info[MAX_NESTING_LEVEL];
__ldm int cur_lev = -1;
static inline int lowbit(int x){
  return x & -x;
}
static void swgomp_level_push(uint64_t bitmap, int tlev){
  struct level_info *info = level_info + (++cur_lev);
  info->bitmap = bitmap;
  info->nthr = __builtin_popcountl(bitmap);
  info->ithr = __builtin_popcountl(bitmap & (1UL << _MYID) - 1);
  info->rmask = bitmap >> (_ROW << 3) & 0xff;
  info->cmask = 0;
  for (int i = 0; i < 8; i ++) {
    if (bitmap & 1UL << ((i << 3) | _COL)) {
      info->cmask |= 1 << i;
    }
  }
  if (bitmap & 0xAAAAAAAAAAAAAAAAUL) {
    info->next_lev = -1;
  } else if (bitmap & 0x4444444444444444UL) {
    info->next_lev = 0;
  } else if (bitmap & 0x1010101010101010UL) {
    info->next_lev = 1;
  } else if (bitmap & 0x0100010001000100UL) {
    info->next_lev = 2;
  } else if (bitmap & 0x0001000000010000UL) {
    info->next_lev = 3;
  } else if (bitmap & 0x0000000100000000UL) {
    info->next_lev = 4;
  } else
    info->next_lev = 5;
  int cpe0 = _MYID;
  while (bitmap >> (cpe0-lowbit(cpe0)) & 1) {
    cpe0 -= lowbit(cpe0);
  }
  info->cpe0 = cpe0;
  info->tlev = tlev;
}
static void swgomp_level_pop(){
  cur_lev --;
}
static uint64_t calc_bitmap(int rootpe, int clev, int nthr) {
  uint64_t bm = 0;
  struct {
    int root;
    int clev;
    int nthr;
  } q[64];
  q[0] = (typeof(q[0])){rootpe, clev, nthr};
  int s = 0, t = 1;
  while (s < t) {
    int root = q[s].root;
    int clev = q[s].clev;
    int nthr = q[s].nthr;
    int nl = (nthr + 1) >> 1;
    int nr = nthr - nl;
    if (nr > 0) {
      q[t++] = (typeof(q[0])){root + (1 << clev), clev - 1, nr};
    }
    if (nl > 1) {
      q[s] = (typeof(q[0])){root, clev - 1, nl};
    } else {
      bm |= (uint64_t)1 << root;
      s++;
    }
  }
  return bm;
}
//#include "slave.h"
static void task_tree_spawn() {
  uint64_t bitmap = tasktree_spawn_arg.s.bitmap;
  int nthr = tasktree_spawn_arg.s.nthr;
  int clev = tasktree_spawn_arg.s.clev;
  int tlev = tasktree_spawn_arg.s.tlev;
  int nl = (nthr + 1) >> 1;
  int nr = nthr - nl;
  int nested = tasktree_spawn_arg.s.nested;
  void (*pc)(void*) = tasktree_spawn_arg.s.pc;
  void *arg = tasktree_spawn_arg.s.arg;
  omp_set_nested(nested);

  if (nr > 0) {
    int child = _MYID + (1 << clev);
    volatile union tree_spawn_arg *child_arg_addr = cpe_pointer(child, &tasktree_spawn_arg);
    union tree_spawn_arg child_arg;
    child_arg.s = (typeof(child_arg.s)){bitmap, clev - 1, tlev, nr, pc, arg, nested};
    child_arg_addr->v = child_arg.v;
    job_spawn_p2p(child, task_tree_spawn, NULL, STK_DO_NOT_MOVE);
  }
  while (nl > 1) {
    clev--;
    nthr = nl;
    nl = (nthr + 1) >> 1;
    nr = nthr - nl;
    int child = _MYID + (1 << clev);
    volatile union tree_spawn_arg *child_arg_addr = cpe_pointer(child, &tasktree_spawn_arg);
    union tree_spawn_arg child_arg;
    child_arg.s = (typeof(child_arg.s)){bitmap, clev - 1, tlev, nr, pc, arg, nested};
    child_arg_addr->v = child_arg.v;
    job_spawn_p2p(child, task_tree_spawn, NULL, STK_DO_NOT_MOVE);
  }
  swgomp_level_push(bitmap, tlev);
  flush_slave_cache();
  task_tree_sync();
  // printf("%p %p\n", pc, arg);
  pc(arg);
  flush_slave_cache();
  task_tree_sync();
  swgomp_level_pop();
}
extern __ldm uint64_t team_bitmap;
void task_tree_start(int nthr, void (*pc)(void *), void *arg){
  if (cur_lev == -1) {
    swgomp_level_push(team_bitmap, 5);
  }
  tasktree_spawn_arg.s.bitmap = calc_bitmap(_MYID, level_info[cur_lev].next_lev, nthr);
  tasktree_spawn_arg.s.nthr = nthr;
  tasktree_spawn_arg.s.clev = level_info[cur_lev].next_lev;
  tasktree_spawn_arg.s.tlev = level_info[cur_lev].next_lev;
  tasktree_spawn_arg.s.pc = pc;
  tasktree_spawn_arg.s.arg = arg;
  tasktree_spawn_arg.s.nested = omp_get_nested();
  long fp;
  task_tree_spawn();
}

void task_tree_sync(){
  __asm__("synr %0\n\t"
          "sync %1\n\t"
          "synr %0\n\t" 
          :: "r"(level_info[cur_lev].rmask), "r"(level_info[cur_lev].cmask)
          : "memory");
}
// 0, 32, 16, 48,  8, 24, 40, 56
// 0,  1,  2,  3,  4,  5,  6,  7
#endif