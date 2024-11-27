#ifndef LWPF_SHARED_CPE_
#define LWPF_SHARED_CPE_
#include <simd.h>
#include <slave.h>

extern __thread_local LWPF_PCVEC lwpf_local_counter[];
#define LWPF_WEAK_DATA __attribute__((weak)) __attribute__((section(".data")))
#define LWPF_WEAK __attribute__((weak))
static void lwpf_sync_counters_m2c(void *lwpf_cpe_counter, int kernel_count) {
  volatile int reply = 0;
  athread_get(PE_MODE, lwpf_cpe_counter, lwpf_local_counter, sizeof(LWPF_PCVEC) * kernel_count, (void*)&reply, 0, 0, 0);
  while (reply != 1);
}

static void lwpf_sync_counters_c2m(void *lwpf_cpe_counter, int kernel_count) {
  //void *local_counter;
  //weak_local(local_counter, "lwpf_local_counter");
  CRTS_dma_put(lwpf_cpe_counter, lwpf_local_counter, sizeof(LWPF_PCVEC) * kernel_count);
  //volatile int reply = 0;
  //athread_put(PE_MODE, lwpf_local_counter, lwpf_cpe_counter, sizeof(LWPF_PCVEC) * kernel_count, (void*)&reply, 0, 0);
  //while (reply != 1);
}

#define lwpf_enter(x) {                                                 \
    lwpf_sync_counters_m2c(lwpf_global_counter_ ## x[_MYID], lwpf_kernel_count_ ## x); \
    evt_conf_t lconf;\
    volatile int reply = 0;\
    athread_get(PE_MODE, &lwpf_evt_config_ ## x, &lconf, sizeof(evt_conf_t), (void*)&reply, 0, 0, 0);\
    while (reply != 1);\
    config_pcrs(lconf.evt);                             \
  }

#define lwpf_exit(x) lwpf_sync_counters_c2m(lwpf_global_counter_ ## x[_MYID], lwpf_kernel_count_ ## x);

#define lwpf_start(kernel) {						\
    LWPF_PCVEC cntrs;							\
    asm volatile(read_pcrs("%1")					\
		 "vsubl %2, %1, %0\n\t"					\
		 : "=" PCVEC(lwpf_local_counter[kernel]), "=&" PCVEC(cntrs) \
		 : "0"(lwpf_local_counter[kernel]));                    \
  }
#define lwpf_stop(kernel) {						\
    LWPF_PCVEC cntrs;                                                   \
    asm volatile(read_pcrs("%1")                                        \
                 "vaddl %2, %1, %0\n\t"                                 \
                 : "=" PCVEC(lwpf_local_counter[kernel]), "=&" PCVEC(cntrs) \
                 : "0"(lwpf_local_counter[kernel]));                    \
  }

#define lwpf_start_volatile(kernel) {					\
    LWPF_PCVEC cntrs;							\
    asm volatile(read_pcrs("%1")					\
		 "vsubl %2, %1, %0\n\t"					\
		 : "=" PCVEC(lwpf_local_counter[kernel]), "=&" PCVEC(cntrs) \
		 : "0"(lwpf_local_counter[kernel]) : "memory");		\
  }
#define lwpf_stop_volatile(kernel) {					\
    LWPF_PCVEC cntrs;                                                   \
    asm volatile_volatile(read_pcrs("%1")				\
			  "vaddl %2, %1, %0\n\t"			\
			  : "=" PCVEC(lwpf_local_counter[kernel]), "=&" PCVEC(cntrs) \
			  : "0"(lwpf_local_counter[kernel]) : "memory"); \
  }

#endif

#ifdef __sw_slave__

#define K(x) x,
typedef enum {
  LWPF_KERNELS
#define U(x) LWPF_KERNELS_END_ ## x
  LWPF_UNIT
#undef U
#define U(x) lwpf_kernel_ ## x
} LWPF_UNIT;
#undef U
#undef K

#define K(x) #x,
#define U(x) lwpf_kernel_names_ ## x
LWPF_WEAK const char *LWPF_UNIT[] = {LWPF_KERNELS "LWPF_KERNELS_END"};
#undef K
#undef U

#define U(x) lwpf_kernel_count_ ## x
extern LWPF_WEAK const long LWPF_UNIT = 
#undef U
#define U(x) LWPF_KERNELS_END_ ## x
  LWPF_UNIT;
#undef U
#define U(x) lwpf_evt_config_ ## x
LWPF_WEAK_DATA evt_conf_t LWPF_UNIT;
#undef U
#define U(x) lwpf_global_counter_ ## x
LWPF_WEAK_DATA long LWPF_UNIT[64][
#undef U
#define U(x) LWPF_KERNELS_END_ ## x
                   LWPF_UNIT
                   ][NPCR];


__thread_local LWPF_PCVEC lwpf_local_counter[LWPF_UNIT] LWPF_WEAK;
#undef U

#define U(x) lwpf_init_ ## x
EXTERNC LWPF_WEAK void LWPF_UNIT(evt_conf_t *conf){
#undef U
  int i;
  evt_conf_t lconf;
  if (!conf){
#define U(x) lwpf_evt_config_ ## x
    lconf.pc_mask = 0;
#ifdef EVT_PC0
    lconf.evt[0] = EVT_PC0;
    lconf.pc_mask |= MASK_PC(0);
#endif
#ifdef EVT_PC1
    lconf.evt[1] = EVT_PC1;
    lconf.pc_mask |= MASK_PC(1);
#endif
#ifdef EVT_PC2
    lconf.evt[2] = EVT_PC2;
    lconf.pc_mask |= MASK_PC(2);
#endif
#ifdef EVT_PC3
    lconf.evt[3] = EVT_PC3;
    lconf.pc_mask |= MASK_PC(3);
#endif
#ifdef EVT_PC4
    lconf.evt[4] = EVT_PC4;
    lconf.pc_mask |= MASK_PC(4);
#endif
#ifdef EVT_PC5
    lconf.evt[5] = EVT_PC5;
    lconf.pc_mask |= MASK_PC(5);
#endif
#ifdef EVT_PC6
    lconf.evt[6] = EVT_PC6;
    lconf.pc_mask |= MASK_PC(6);
#endif
#ifdef EVT_PC7
    lconf.evt[7] = EVT_PC7;
    lconf.pc_mask |= MASK_PC(7);
#endif
  } else {
    for (i = 0; i < NPCR; i ++)
      lconf.evt[i] = conf->evt[i];
    lconf.pc_mask = conf->pc_mask;
  }
  if (!_MYID) {
    volatile int reply = 0;
    athread_put(PE_MODE, &lconf, &(LWPF_UNIT), sizeof(evt_conf_t), (void*)&reply, 0, 0);
    while (reply != 1);
  }
//if (conf == NULL) puts("conf == NULL");
//else printf("conf = %x\n", conf);
//
//printf("wyj_cpe mask = 0x%x\n", LWPF_UNIT.pc_mask);
//for (i = 0; i < 8; i++)
//    printf("MYID= %d, LWPF_UNIT.evt[%d] = 0x%x\n", _MYID, i, LWPF_UNIT.evt[i]);


#undef U
#define U(x) LWPF_KERNELS_END_ ## x
  LWPF_PCVEC v0 = 0;

  for (i = 0; i < LWPF_UNIT; i ++){
    lwpf_local_counter[i] = v0;
  }
#undef U
#define U(x) lwpf_exit(x);
  LWPF_UNIT
#undef U
}
#endif

