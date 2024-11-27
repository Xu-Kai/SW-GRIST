#pragma once
#include "swgomp-stack.h"
void jobserver();
void job_spawn_p2p(int tid, void (*pc)(void*), void *arg, enum stack_loc stkloc);
void job_spawn_bcast(int mode, int mask, void (*pc)(void*), void *arg, enum stack_loc stkloc) ;
void job_wait_p2p(int tid);
void slave_jobserver();
void jobserver_athread_init();
void jobserver_athread_spawn();