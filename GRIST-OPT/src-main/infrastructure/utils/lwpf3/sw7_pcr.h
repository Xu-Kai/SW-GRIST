#ifndef SW7_PCR_H_
#define SW7_PCR_H_
#include <stdlib.h>
#include <string.h>

#define PCR_DEF_BEGIN(name) enum name ## _OPTS {
#define PCR_OPT_DEF(opt, value) opt = value,
#define PCR_DEF_END(name) name ## _MAX};
#include "sw7_pcrdef.h"
#undef PCR_DEF_BEGIN
#undef PCR_OPT_DEF
#undef PCR_DEF_END

#define PCR_DEF_BEGIN(name) static char *name ##_NAMES[] = {
#define PCR_OPT_DEF(opt, value) #opt,
#define PCR_DEF_END(name) #name "_END"};
#include "sw7_pcrdef.h"
#undef PCR_DEF_BEGIN
#undef PCR_OPT_DEF
#undef PCR_DEF_END
#define PCR0 0x20
#define PCR1 0x21
#define PCR2 0x22
#define PCR3 0x23
#define PCR4 0x24
#define PCR5 0x25
#define PCR6 0x26
#define PCR7 0x27

#define PCRC 0x8
#define PCRC_ALL 0xff
#define NPCR 8
static char **PCR_NAMES[] = {PC0_NAMES, PC1_NAMES, PC2_NAMES, PC3_NAMES,
		      PC4_NAMES, PC5_NAMES, PC6_NAMES, PC7_NAMES};
#define PCR_SEL_OFF 56

#define PC_ALL 0xff
#define read_pcr_internal(target, num) "rcsr " target ", " #num "\n\t"
#define read_pcr(target, num) read_pcr_internal(target, num)

#define read_pcr0(target) read_pcr(target, PCR0)
#define read_pcr1(target) read_pcr(target, PCR1)
#define read_pcr2(target) read_pcr(target, PCR2)
#define read_pcr3(target) read_pcr(target, PCR3)
#define read_pcr4(target) read_pcr(target, PCR4)
#define read_pcr5(target) read_pcr(target, PCR5)
#define read_pcr6(target) read_pcr(target, PCR6)
#define read_pcr7(target) read_pcr(target, PCR7)

#define insf0(x, y) "vinsf " x ", " y ", 0, " y "\n\t"
#define insf1(x, y) "vinsf " x ", " y ", 1, " y "\n\t"
#define insf2(x, y) "vinsf " x ", " y ", 2, " y "\n\t"
#define insf3(x, y) "vinsf " x ", " y ", 3, " y "\n\t"
#define insf4(x, y) "vinsf " x ", " y ", 4, " y "\n\t"
#define insf5(x, y) "vinsf " x ", " y ", 5, " y "\n\t"
#define insf6(x, y) "vinsf " x ", " y ", 6, " y "\n\t"
#define insf7(x, y) "vinsf " x ", " y ", 7, " y "\n\t"

#ifdef __sw_slave__
static inline void config_pcrs(long *conf){
  asm("wcsr %0, %1\n\t" :: "r"(PCRC_ALL), "i"(PCRC));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[0]) << PCR_SEL_OFF), "i"(PCR0));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[1]) << PCR_SEL_OFF), "i"(PCR1));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[2]) << PCR_SEL_OFF), "i"(PCR2));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[3]) << PCR_SEL_OFF), "i"(PCR3));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[4]) << PCR_SEL_OFF), "i"(PCR4));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[5]) << PCR_SEL_OFF), "i"(PCR5));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[6]) << PCR_SEL_OFF), "i"(PCR6));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[7]) << PCR_SEL_OFF), "i"(PCR7));
}

#define read_pcrs(cntrs)                                \
  "rcsr " cntrs ", 0x27\n\t"				\
  "vinsf " cntrs ", " cntrs ", 7, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x26\n\t"				\
  "vinsf " cntrs ", " cntrs ", 6, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x25\n\t"				\
  "vinsf " cntrs ", " cntrs ", 5, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x24\n\t"				\
  "vinsf " cntrs ", " cntrs ", 4, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x23\n\t"				\
  "vinsf " cntrs ", " cntrs ", 3, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x22\n\t"				\
  "vinsf " cntrs ", " cntrs ", 2, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x21\n\t"				\
  "vinsf " cntrs ", " cntrs ", 1, " cntrs "\n\t"        \
  "rcsr " cntrs ", 0x20\n\t"				\

#define LWPF_PCVEC int512
#define PCVEC "f"
#endif
#endif
