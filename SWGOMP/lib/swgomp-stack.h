#pragma once
#ifndef __ASSEMBLER__
enum stack_loc {
  STK_IN_LDM,
  STK_IN_SHR,
  STK_IN_PRIV,
  STK_DO_NOT_MOVE
};

void shared_stack_init();
void *shared_stack_get(int pe);
void set_stack_run(void (*pc)(void*), void *arg, void *stk, void **oldstk);
void switch_stack_run(void (*pc)(void*), void *arg, enum stack_loc loc);

#ifdef __sw_slave__
void init_saved_sp();
#endif
extern void *swgomp_shared_stack;
#else
#define STK_IN_LDM      0
#define STK_IN_SHR      1
#define STK_IN_PRIV     2
#define STK_DO_NOT_MOVE 3
#endif