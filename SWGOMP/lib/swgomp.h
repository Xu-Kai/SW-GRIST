#pragma once
int swgomp_get_cpe_num();
int swgomp_get_rootpe();
void GOMP_barrier();
void GOMP_parallel(void (*fn)(void*), void *arg, unsigned nthr, unsigned flags);