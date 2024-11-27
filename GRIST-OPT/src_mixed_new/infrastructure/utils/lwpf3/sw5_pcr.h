#ifndef SW5_PCR_H_
#define SW5_PCR_H_
#include <stdlib.h>
#include <string.h>

#define PCR_DEF_BEGIN(name) enum name ## _OPTS {
#define PCR_OPT_DEF(opt, value) opt = value,
#define PCR_DEF_END(name) name ## _MAX};
#include "sw5_pcrdef.h"
#undef PCR_DEF_BEGIN
#undef PCR_OPT_DEF
#undef PCR_DEF_END

#define PCR_DEF_BEGIN(name) static char *name ##_NAMES[] = {
#define PCR_OPT_DEF(opt, value) #opt,
#define PCR_DEF_END(name) #name "_END"};
#include "sw5_pcrdef.h"
#undef PCR_DEF_BEGIN
#undef PCR_OPT_DEF
#undef PCR_DEF_END

#define PCR0 0x5
#define PCR1 0x6
#define PCR2 0x7
#define PCR3 0x4
#define PCRC 0x8

#define PCRC_ALL_FLOP          (0x3f << 16)
#define PCRC_ALL_PC            (0x1f)
#define PCRC_ALL               (PCRC_ALL_FLOP | PCRC_ALL_PC)
#define PCR_SEL_OFF            59
#define NPCR 4
static char **PCR_NAMES[] = {PC0_NAMES, PC1_NAMES, PC2_NAMES, PC3_NAMES};
//#define PCRC_VDIV_VSQRT        (1 << 21) //浮点向量除法平方根指令计数使能。本位为1，从核性能计数器1的计数使能为1，且计数事件为5’h15（浮点除法、平方根类指令等效操作计数）时，浮点向量除法平方指令参与计数（每次+4）。
//#define PCRC_FDIV_FSQRT        (1 << 20) //浮点标量除法平方根指令计数使能。本位为1，从核性能计数器1的计数使能为1，且计数事件为5’h15（浮点除法、平方根类指令等效操作计数）时，浮点标量除法平方指令参与计数（每次+1）。
//#define PCRC_VMA               (1 << 19) //浮点向量乘加类指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点向量乘加类指令参与计数（每次+8）。
//#define PCRC_VADD_VSUB_VMUL    (1 << 18) //浮点向量加减乘指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点向量加减乘指令参与计数（每次+4）。
//#define PCRC_FMA               (1 << 17) //浮点标量乘加类指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点标量乘加类指令参与计数（每次+2）。
//#define PCRC_FADD_FSUB_FMUL    (1 << 16) //浮点标量加减乘指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点标量加减乘指令参与计数（每次+1）。
#define PCRC_PC2_OVERFLOW      (1 <<  5) //TA性能计数器的溢出中断使能。为1时使能。   
#define PCRC_PC2               (1 <<  4) //TA性能计数器的计数使能。为1时使能。        
#define PCRC_PC1_OVERFLOW      (1 <<  3) //从核性能计数器1的溢出中断使能。为1时使能。
#define PCRC_PC1               (1 <<  2) //从核性能计数器1的计数使能。为1时使能。    
#define PCRC_PC0_OVERFLOW      (1 <<  1) //从核性能计数器0的溢出中断使能。为1时使能。
#define PCRC_PC0               (1 <<  0) //从核性能计数器0的计数使能。为1时使能。


#ifdef __sw_slave__ //The following macros should be only available on CPEs!!!
static inline void config_pcrs(long *conf){
  asm("wcsr %0, %1\n\t" :: "r"(PCRC_ALL), "i"(PCRC));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[0]) << 59), "i"(PCR0));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[1]) << 59), "i"(PCR1));
  asm("wcsr %0, %1\n\t" :: "r"(((long)conf[2]) << 59), "i"(PCR2));
}

#define read_pcrs(cntrs)                                \
  "rcsr " cntrs ", 4\n\t"                               \
  "vinsf " cntrs ", " cntrs ", 3, " cntrs "\n\t"        \
  "rcsr " cntrs ", 7\n\t"                               \
  "vinsf " cntrs ", " cntrs ", 2, " cntrs "\n\t"        \
  "rcsr " cntrs ", 6\n\t"                               \
  "vinsf " cntrs ", " cntrs ", 1, " cntrs "\n\t"        \
  "rcsr " cntrs ", 5\n\t"

#define LWPF_PCVEC int256
#define PCVEC "r"
#endif
#endif

