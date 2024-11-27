#pragma once
#include "swdefs.h"
struct level_info {
  uint64_t bitmap;
  int tlev;
  int next_lev;
  int nthr;
  int ithr;
  int cpe0;
  char rmask, cmask;
};
extern __ldm int cur_lev;
extern __ldm struct level_info level_info[6];
extern void task_tree_start(int nthr, void (*pc)(void*), void *arg);
extern void task_tree_sync();