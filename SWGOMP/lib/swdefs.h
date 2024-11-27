#pragma once
#include <stdint.h>
#include <stddef.h>
extern int cal_locked_printf(const char *, ...);

#define MAX_NESTING_LEVEL 6
#define MAX_NUM_PES 64
typedef int int512 __attribute__((__mode__(__V16SI__)));
#ifdef __sw_slave__
typedef int desc_t __attribute__((__mode__(__V16SI__)));
#endif
#define __ldm __attribute__((section(".ldm")))
extern __ldm int _MYID;
extern __ldm char _ROW, _COL;

#ifndef RMA_MODE_SINGLE_CORE
#define RMA_MODE_SINGLE_CORE                    0ULL
#define RMA_MODE_ROW_BCAST                      1ULL
#define RMA_MODE_COL_BCAST                      2ULL
#define RMA_MODE_ROW_BCAST_LOCAL_RECEIVE_DATA   5ULL
#define RMA_MODE_COL_BCAST_LOCAL_RECEIVE_DATA   6ULL

#define RMA_OPCODE_PUT          0ULL
#define RMA_OPCODE_GET          1ULL
#define RMA_OPCODE_BARRIER      5ULL
#define RMA_OPCODE_ALL_BARRIER  6ULL
#endif

static inline void unsafe_wait_value(volatile long *rply, long target) {
  long rplyv;
  __asm__("1:\n\t"
          "ldl %0, %1\n\t"
          "subl %0,%2,%0\n\t"
          "bne %0, 1b\n\t"
          "memb\n\t"
          :"=&r"(rplyv), "+m"(*rply) : "r"(target) : "memory");
}

#ifdef __sw_slave__
extern void flush_slave_cache();
static inline long faal(long *ptr) {
  long ret;
  asm volatile("faal %0, 0(%2)\n\t" : "=r"(ret), "+m"(*ptr):"r"(ptr));
  return ret;
}
static inline void rma(int R_tid, void *Laddr, void *Raddr, size_t size, volatile long *Lrply, volatile long *Rrply, int mode, int opcode, int mask){
  uint64_t Rb, Rc;
  union {
    uint64_t l[2];
    int512 v;
  } desc;
  desc.l[1]=(uint64_t)(((((uint64_t)(Rrply))&0x3ffffUL)<<32)
                    |	(((uint64_t)(mask))<<24)
                    | (((uint64_t)(Lrply))&0x3ffff));
  desc.l[0]=(uint64_t)((((uint64_t)(mode))<<60)
                    | (((uint64_t)(opcode))<<56) | (uint64_t)(size));
  Rb = (((uint64_t)R_tid)<<20)
              | (((uint64_t)(Raddr))&0x3ffffUL);
  Rc = ((uint64_t)(Laddr))&0x3ffffUL;
  __asm__ volatile("rma %0,%1,%2\n\tmemb\n\t" :"+r"(desc.v),"+r"(Rb),"+r"(Rc),"+m"(*Lrply)::"memory");
}
static __ldm volatile long lrply = 0;
static inline void rma_bcast(void *addr, size_t size, volatile long *rply, int mode, int mask){
  lrply = 0;
  // cal_locked_printf("bcast: %p %d %p %d %d\n", addr, size, rply, mode, mask);
  rma(0, addr, addr, size, &lrply, rply, mode, RMA_OPCODE_PUT, mask);
  unsafe_wait_value(&lrply, 1);
}

static inline void *cpe_pointer(int tid, void *ptr){
  uint64_t remote_addr = (uint64_t)ptr | 1UL << 45 | (long)tid << 20;
  return (void*)remote_addr;
}
static inline void *cpe_symbol_addr(void *arg) {
  return arg;
}
#else
static inline int arrayid(void)
{
    unsigned long cid; 
    asm volatile(
            "rcid %0\n"
            "sll  %0, 61, %0\n"
            "srl  %0, 61, %0\n"
            :"=r"(cid));
    return (int) cid;
}
static inline void *ldm_addr(int cgid, int speid, void *disp){
  uint64_t ldm_hw_base = 0x1ULL<<47ULL | 0x2ULL<<22ULL;
  uint64_t addr = (uint64_t)disp;
  return (void*)(ldm_hw_base + ( (unsigned long)cgid<<40) + ((unsigned long)speid <<24) + addr);
}

static inline void *cpe_pointer(int tid, void *ptr){
  // printf("%p\n", ptr);
  return ldm_addr(arrayid(), tid, ptr);
  // return (void*)remote_addr;
}
static const uint64_t ldm_offset = 0x500000004000UL;
static inline void *cpe_symbol_addr(void *arg) {
  // printf("cpe_symbol: %p\n", arg);
  return (void*)((uint64_t)arg - ldm_offset);
}
#endif