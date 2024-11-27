#define NCO 128
#define CO_PE 8
#define CO_PP (NCO/CO_PE)
#define NCI 128
#define L 32
#define NBATCH 16
#define BA_PE (64/CO_PE)
#define VWID 8
#define CIBLK 8
#include <math.h>
struct conv1d3_128x128x32_param {
  float (*conv_out)[NCO][L];
  float (*conv_in)[NCI][L];
  float (*kern)[NCI][3][VWID];
  float *bias;
  int nbatch;
};
#ifdef __sw_slave__
#define LWPF_KERNELS K(RA) K(RK) K(RB) K(W) K(M) K(A)
#define LWPF_UNIT U(CONV1D)
#define EVT_PC0 PC0_CYCLE
#define EVT_PC1 PC1_INST
#define EVT_PC3 PC3_L1IC_MISS
#define EVT_PC4 PC4_L1IC_MISSTIME
#define EVT_PC6 PC6_INST_VECTOR_FLOAT_ADDS
#define EVT_PC7 PC7_INST_VECTOR_FLOAT_MULS
#include "lwpf3/lwpf.h"
#include <simd.h>
#include <slave.h>
void myprintdv8(const char *name, doublev8 a){
  printf(name);
  double *b = &a;
  for (int i = 0; i < 8; i ++) {
    printf("%g ", b[i]);
  }
  puts("");
}
#define printdv8(x) myprintdv8(#x, x);
void myprintfv8(const char *name, floatv8 a){
  printf(name);
  float *b = &a;
  for (int i = 0; i < 8; i ++) {
    printf("%g ", b[i]);
  }
  printf("\n");
}
#define printfv8(x) myprintfv8(#x": ", x);
__always_inline void dma_iget(void *mem, void *ldm, size_t count, volatile long *rply) {
  //long dma_desc[4];
  int512 dma_desc;
  asm volatile(
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "ldi %[DESC], 1($31)\n\t"
      "sll %[DESC], 56, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      : [ DESC ] "=&r"(dma_desc), "+m"(*rply)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(rply), [ SIZE ] "r"(count)
      : "memory");
}

__always_inline void dma_iput(void *mem, void *ldm, size_t count, volatile long *rply) {
  int512 dma_desc;
  asm volatile(
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc), "+m"(*rply)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(rply), [ SIZE ] "r"(count)
      : "memory");
}

__always_inline void unsafe_wait_value(volatile long *rply, long target) {
  long rplyv;
  __asm__ volatile("1:\n\t"
                   "ldl %0, %1\n\t"
                   "subl %0,%2,%0\n\t"
                   "bne %0, 1b\n\t"
                   "memb\n\t"
                   :"=&r"(rplyv), "+m"(*rply) : "r"(target) : "memory");
}
__always_inline void unsafe_wait_value_ge(volatile long *rply, long target) {
  long rplyv;
  __asm__ volatile("1:\n\t"
                   "ldl %0, %1\n\t"
                   "subl %0,%2,%0\n\t"
                   "blt %0, 1b\n\t"
                   "memb\n\t"
                   :"=&r"(rplyv), "+m"(*rply) : "r"(target) : "memory");
}
__asm__(
  ".macro maybe  cond, inst:vararg\n\t"
  "  .if \\cond\n\t"
  "    \\inst\n\t"
  "  .endif\n\t"
  ".endm\n\t"

  ".macro iter1 curr, stlast=1, ldnext=1, out=$0, in=$1, p8=$2, m8=$3, nxtr=$34\n\t"
  "  vlenma $36, $32, $40, $40 ; maybe \\stlast, vsts $50,   0+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $41, $41 ; maybe \\stlast, vsts $51, 128+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $42, $42 ; maybe \\stlast, vsts $52, 256+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $43, $43 ; maybe \\stlast, vsts $53, 384+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $44, $44 ; maybe \\stlast, vsts $54, 512+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $45, $45 ; maybe \\stlast, vsts $55, 640+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $46, $46 ; maybe \\stlast, vsts $56, 768+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $47, $47 ; maybe \\stlast, vsts $57, 896+\\curr-32(\\out)\n\t"

  "  vlenma $35, $38, $40, $40 ; maybe \\ldnext, vlds $50,   0+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $41, $41 ; maybe \\ldnext, vlds $51, 128+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $42, $42 ; maybe \\ldnext, vlds $52, 256+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $43, $43 ; maybe \\ldnext, vlds $53, 384+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $44, $44 ; maybe \\ldnext, vlds $54, 512+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $45, $45 ; maybe \\ldnext, vlds $55, 640+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $46, $46 ; maybe \\ldnext, vlds $56, 768+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $47, $47 ; maybe \\ldnext, vlds $57, 896+\\curr+32(\\out)\n\t"

  "  vlenma $37, $39, $40, $40 ; maybe \\ldnext, vlds $33, \\curr(\\in)\n\t"
  "  vlenma $37, $39, $41, $41 ; maybe \\ldnext, vlds $32, \\curr+32(\\in)\n\t"
  "  vlenma $37, $39, $42, $42 ; maybe \\ldnext, vlds \\nxtr, \\curr+64(\\in)\n\t"
  "  vlenma $37, $39, $43, $43 \n\t"
  "  vlenma $37, $39, $44, $44 \n\t"
  "  vlenma $37, $39, $45, $45 \n\t"
  "  vlenma $37, $39, $46, $46 ; maybe \\ldnext vconw  $32, $33, \\m8, $38\n\t"
  "  vlenma $37, $39, $47, $47 ; maybe \\ldnext vconw \\nxtr, $32, \\p8, $39\n\t"
  ".endm\n\t"

  ".macro iter2 curr, stlast=1, ldnext=1, out=$0, in=$1, p8=$2, m8=$3, nxtr=$34\n\t"
  "  vlenma $36, $32, $50, $50 ; maybe \\stlast, vsts $40,   0+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $51, $51 ; maybe \\stlast, vsts $41, 128+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $52, $52 ; maybe \\stlast, vsts $42, 256+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $53, $53 ; maybe \\stlast, vsts $43, 384+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $54, $54 ; maybe \\stlast, vsts $44, 512+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $55, $55 ; maybe \\stlast, vsts $45, 640+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $56, $56 ; maybe \\stlast, vsts $46, 768+\\curr-32(\\out)\n\t"
  "  vlenma $36, $32, $57, $57 ; maybe \\stlast, vsts $47, 896+\\curr-32(\\out)\n\t"

  "  vlenma $35, $38, $50, $50 ; maybe \\ldnext, vlds $40,   0+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $51, $51 ; maybe \\ldnext, vlds $41, 128+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $52, $52 ; maybe \\ldnext, vlds $42, 256+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $53, $53 ; maybe \\ldnext, vlds $43, 384+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $54, $54 ; maybe \\ldnext, vlds $44, 512+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $55, $55 ; maybe \\ldnext, vlds $45, 640+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $56, $56 ; maybe \\ldnext, vlds $46, 768+\\curr+32(\\out)\n\t"
  "  vlenma $35, $38, $57, $57 ; maybe \\ldnext, vlds $47, 896+\\curr+32(\\out)\n\t"

  "  vlenma $37, $39, $50, $50 ; maybe \\ldnext, vlds $33, \\curr(\\in)\n\t"
  "  vlenma $37, $39, $51, $51 ; maybe \\ldnext, vlds $32, \\curr+32(\\in)\n\t"
  "  vlenma $37, $39, $52, $52 ; maybe \\ldnext, vlds \\nxtr, \\curr+64(\\in)\n\t"
  "  vlenma $37, $39, $53, $53 \n\t"
  "  vlenma $37, $39, $54, $54 \n\t"
  "  vlenma $37, $39, $55, $55 \n\t"
  "  vlenma $37, $39, $56, $56 ; maybe \\ldnext vconw  $32, $33, \\m8, $38\n\t"
  "  vlenma $37, $39, $57, $57 ; maybe \\ldnext vconw \\nxtr, $32, \\p8, $39\n\t"
  ".endm\n\t"
);
__always_inline void conv1d3x8x32(float *out, float *in, float *kern){
  int i;
  __asm__ volatile(
  "vlds   $32, 0(%[IN])   ; "
  "vlds   $34, 32(%[IN])  ; vcpyf  $31, $33\n\t"
  "vlds   $35, 0(%[KERN])     ; \n\t"
  "vlds   $36, 32(%[KERN])    ; "
  "vlds   $37, 64(%[KERN])    ; wcsr   %[T0], 0x92\n\t"
  "vconw  $32, $33, %[M8], $38\n\t"
  "vlds   $40,    0(%[OUT])\n\t"
  "vconw  $34, $32, %[P8], $39\n\t"

  "vlds $41, 128(%[OUT])\n\t"
  "vlds $42, 256(%[OUT])\n\t"
  "vlds $43, 384(%[OUT])\n\t"
  "vlds $44, 512(%[OUT])\n\t"
  "vlds $45, 640(%[OUT])\n\t"
  "vlds $46, 768(%[OUT])\n\t"
  "vlds $47, 896(%[OUT])\n\t"

  "iter1   0, 0, 1, %[OUT], %[IN], %[P8], %[M8], $34\n\t"
  "iter2  32, 1, 1, %[OUT], %[IN], %[P8], %[M8], $34\n\t"
  "iter1  64, 1, 1, %[OUT], %[IN], %[P8], %[M8], $34\n\t"
  "iter2  96, 1, 0, %[OUT], %[IN], %[P8], %[M8], $31\n\t"
  "1:\n\t"

  "vsts $50,   0+96(%[OUT])\n\t"
  "vsts $51, 128+96(%[OUT])\n\t"
  "vsts $52, 256+96(%[OUT])\n\t"
  "vsts $53, 384+96(%[OUT])\n\t"
  "vsts $54, 512+96(%[OUT])\n\t"
  "vsts $55, 640+96(%[OUT])\n\t"
  "vsts $56, 768+96(%[OUT])\n\t"
  "vsts $57, 896+96(%[OUT])\n\t"
  :[OUT]"+&r"(out), [IN]"+&r"(in), [KERN]"+&r"(kern), [I]"=&r"(i): [P8]"r"(8L), [M8]"r"(-8L), [T0]"r"(8L << 32)
  : "$32", "$33", "$34", "$35", "$36", "$37", "$38", "$39", "$40", 
    "$41", "$42", "$43", "$44", "$45", "$46", "$47", "$50", "$51",
    "$52", "$53", "$54", "$55", "$56", "$57", "memory"
  );
}

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

// 
static __ldm_kernel(conv1d) __attribute__((aligned(64))) float conv_in[NCI][L], conv_out[CO_PP][L], kernpp[CO_PP/VWID][NCI][3][VWID];
static __ldm_kernel(conv1d) volatile long s_rma_rply, r_rma_rply, dma_rply, rply, w_rply;
const int x = sizeof(conv_in) + sizeof(conv_out) + sizeof(kernpp);

void conv1d3_128x128(struct conv1d3_128x128x32_param *pm){
  lwpf_enter(CONV1D);
  long dbgs[8];
  int baid = _MYID / BA_PE;
  int coid = _MYID % BA_PE;
  lwpf_start(A);
  rply = 0;
  w_rply = 0;
  // printf("%p %p %d\n", pm->kern[coid*CO_PP], kernpp);
  dma_iget(pm->kern[coid*CO_PP/VWID], kernpp, CO_PP*NCI*3*sizeof(float), &rply);
  unsafe_wait_value(&rply, 1);
  for (int ib = baid; ib < pm->nbatch; ib += BA_PE) {
    lwpf_start(RB);
    r_rma_rply = 0;
    s_rma_rply = 0;
    rply = 0;
    asm volatile("synr %0\n\t" :: "r"(0xff) : "memory");
    if (coid == 0) {
      dma_iget(pm->conv_in[ib][0], conv_in[0], CIBLK*L*sizeof(float), &rply);
    }
    for (int i = 0; i < CO_PP; i ++) {
      int gco = CO_PP*coid + i;
      floatv8 bias = pm->bias[gco];
      for (int j = 0; j < L; j += 8) {
        simd_store(bias, &conv_out[i][j]);
      }
    }
    // if (_MYID == 0 && ib == 0) {
    //   for (int i = 0; i < 16; i ++){
    //     printf("%g ", conv_out[0][i]);
    //   }
    //   puts("");
    // }
    lwpf_stop(RB);
    w_rply = 0;
    for (int Ci = 0; Ci < NCI; Ci += CIBLK) {
      if (coid == 0) {
        lwpf_start(RA);
        unsafe_wait_value(&rply, Ci/CIBLK+1);
        lwpf_stop(RA);

        rma(0, conv_in[Ci], conv_in[Ci], CIBLK*L*sizeof(float), &s_rma_rply, &r_rma_rply, RMA_MODE_ROW_BCAST, RMA_OPCODE_PUT, 1<<_ROW);
        unsafe_wait_value(&s_rma_rply, Ci/CIBLK+1);
        asm volatile("synr %0\n\t" :: "r"(0xff) : "memory");
        // asm volatile("synr %0\n\t" :: "r"(0xff) : "memory");
        if (Ci < NCI-CIBLK) dma_iget(pm->conv_in[ib][Ci+CIBLK], conv_in[Ci+CIBLK], CIBLK*L*sizeof(float), &rply);
      } else {
        lwpf_start(RK);
        asm volatile("synr %0\n\t" :: "r"(0xff) : "memory");
        unsafe_wait_value(&r_rma_rply, Ci/CIBLK+1);
        lwpf_stop(RK);
      }
      lwpf_start(M);

      for (int Co = 0; Co < CO_PP; Co += VWID) {
        int gco = coid * CO_PP;
        for (int Cii = 0; Cii < CIBLK; Cii ++) {
          conv1d3x8x32(conv_out[Co], conv_in[Ci+Cii], kernpp[Co/VWID][Ci+Cii][0]);
        }
        if (Ci == NCI - CIBLK) {
          for (int i = 0; i < VWID; i ++)
            conv_out[Co+i][30] = conv_out[Co+i][31] = 0;
          dma_iput(pm->conv_out[ib][coid*CO_PP+Co], conv_out[Co], VWID*L*sizeof(float), &w_rply);
        }
      }
      lwpf_stop(M);
    }

    // w_rply = 0;
    lwpf_start(W);
    unsafe_wait_value(&w_rply, CO_PP / VWID);
    // unsafe_wait_value(&w_rply, CO_PP / VWID);
    lwpf_stop(W);
  }
  lwpf_stop(A);
  lwpf_exit(CONV1D);
}
#pragma omp declare target(conv1d3_128x128)

#endif
#ifdef __sw_host__
#include <athread.h>
#define LWPF_UNITS U(CONV1D)
#include "lwpf3/lwpf.h"
extern void slave_conv1d3_128x128();

__attribute__((aligned(256))) float kerno[NCI][NCO][3], bias[NCO], ref[NBATCH][NCO][L];
__attribute__((aligned(256))) float ga[NBATCH][NCI][L], gb[NBATCH][NCO][L], gkern[NCO/VWID][NCI][3][VWID];
void convert_weight(float (*out)[NCI][3][VWID], float (*in)[NCI][3]) {
  for (int i = 0; i < NCO; i ++) {
    for (int j = 0; j < NCI; j ++) {
      for (int k = 0; k < 3; k ++) {
        out[i/VWID][j][k][i%VWID] = -in[i][j][k];
      }
    }
  }
}
#ifdef CONV1D_TEST
int main(){
  athread_init();
  FILE *f = fopen("conv1d.bin", "rb");
  fread(kerno, 4, NCI*NCO*3, f);
  fread(bias, 4, NCO, f);
  fread(ga, 4, NBATCH*NCI*L, f);
  fread(ref, 4, NBATCH*NCO*L, f);
  // for (int i = 0; i < NCO; i ++) {
  //   for (int j = 0; j < NCI; j ++) {
  //     for (int k = 0; k < 3; k ++) {
  //       gkern[i/VWID][j][k][i%VWID] = -kerno[i][j][k];
  //     }
  //   }
  // }
  convert_weight(gkern, kerno);
  struct conv1d3_128x128x32_param pm = {gb, ga, gkern, bias, NBATCH};
  lwpf_init(NULL);
  __real_athread_spawn(slave_conv1d3_128x128, &pm, 1);
  athread_join();

  float maxabs = 0;
  for (int ib = 0; ib < 16; ib ++){
    for (int co = 0; co < 128; co ++) {
      for (int iv = 0; iv < 30; iv ++) {
        if (fabs(ref[ib][co][iv] - gb[ib][co][iv]) > maxabs) {
          maxabs = fabs(ref[ib][co][iv] - gb[ib][co][iv]);
          printf("%d %d %d %f %f %f\n", ib, co, iv, maxabs, ref[ib][co][iv], gb[ib][co][iv]);
        }
      }
    }
  }
  printf("abs: %g\n", maxabs);
  float maxrel = 0;
  for (int ib = 0; ib < 16; ib ++){
    for (int co = 0; co < 128; co ++) {
      for (int iv = 0; iv < 30; iv ++) {
        float diff = fabs(ref[ib][co][iv] - gb[ib][co][iv]);
        float max = fabs(ref[ib][co][iv]) > (gb[ib][co][iv]) ? fabs(ref[ib][co][iv]) : (gb[ib][co][iv]);
        if (diff / max > maxrel) {
          maxrel = diff/max;
          printf("%d %d %d %f %f %f\n", ib, co, iv, maxrel, ref[ib][co][iv], gb[ib][co][iv]);
        }
      }
    }
  }
  printf("rel: %g\n", maxrel);
  // puts("");
  // for (int i = 0; i < 16; i ++){
  //   printf("%g ", gb[0][0][i]);
  // }
  // puts("");
  lwpf_report_summary(stdout);

}
#endif
#endif