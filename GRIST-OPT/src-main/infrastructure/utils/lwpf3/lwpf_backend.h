#ifndef LWPF_BE_
#define LWPF_BE_

#if !defined(SW5) && !defined(SW7) && !defined(SW9)
#error "No architecture selected, try #define SW5/SW7/SW9 before including LWPF headers"
#endif
#ifdef SW5
#include "sw5_pcr.h"
#endif

#ifdef SW7
#include "sw7_pcr.h"
#endif

#ifdef SW9
#include "sw9_pcr.h"
#endif
typedef struct evt_conf {
  long pc_mask;
  long evt[NPCR];
} evt_conf_t;

#define MASK_PC(i) (1 << (i))
#endif
