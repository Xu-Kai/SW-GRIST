#pragma once
#include "swdefs.h"
void taskarray_arg_set(int tid, void (*pc)(void*), void *arg, int nthr, int ldmstack);
void taskarray_spawn_teams();
void task_array_spawn_threads(int nthr, void (*fn)(void*), void *dat);
void task_array_sync();
void taskarray_barrier();
struct taskarray_info{
  int cpe0, row_mask;
  int nthr, ithr;
  int nteam, iteam;
};
extern __ldm struct taskarray_info taskarray_info;