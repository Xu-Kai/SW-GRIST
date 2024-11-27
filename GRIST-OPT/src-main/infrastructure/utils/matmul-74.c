#define GB 512
#define GM 512
#define GN 512
#define B 16
#define M 8
#define N 512
#define INSTRIDE "64"
#define WTSTRIDE "32"
#include <limits.h>
#ifdef __sw_slave__
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
#define LWPF_KERNELS K(REF) K(SIMD) K(ALL) K(GET1) K(GET2) K(PUT)
#define LWPF_UNIT U(MM)
#define EVT_PC0 PC0_CYCLE
#define EVT_PC1 PC1_INST
#define EVT_PC3 PC3_L1IC_MISS
#define EVT_PC4 PC4_L1IC_MISSTIME
#define EVT_PC6 PC6_INST_VECTOR_FLOAT_ADDS
#define EVT_PC7 PC7_INST_VECTOR_FLOAT_MULS
#include "lwpf3/lwpf.h"

void debug_csr(){
  long csr90, csr91, csr92;
  asm volatile("rcsr %0, 0x90\n\trcsr %1, 0x91\n\trcsr %2, 0x92\n\t" : "=r"(csr90), "=r"(csr91), "=r"(csr92));
  printf("90:%016lx 91:%016lx 92:%016lx\n", csr90, csr91, csr92);
}
#define DO(x) x
#define DONOT(x)
#define MM_ITER(IN, WT, INNEW, WTNEW, INPTR, WTPTR, LDNEXT)                    \
  "vlenma " WT ", " IN ", $40, $40\n\t" LDNEXT("vlds "INNEW", " INPTR"\n\t") \
  "vlenma " WT ", " IN ", $41, $41\n\t" LDNEXT("vlds "WTNEW", " WTPTR"\n\t") \
  "vlenma " WT ", " IN ", $42, $42\n\t" \
  "vlenma " WT ", " IN ", $43, $43\n\t" \
  "vlenma " WT ", " IN ", $44, $44\n\t" \
  "vlenma " WT ", " IN ", $45, $45\n\t" \
  "vlenma " WT ", " IN ", $46, $46\n\t" \
  "vlenma " WT ", " IN ", $47, $47\n\t" 

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
#endif
struct mm_par_t {
  float (*gout)[GM][B], (*gin)[GN][B], (*gw)[GN][M], *gbias;
  int nB;
};
#ifdef __sw_slave__

void plain_matmul(float (*lout)[B], float (*lin)[B], float (*lnw_tran)[M], float *bias){
  asm volatile("wcsr %0, 0x92" :: "r"(8L<<32) : "memory");
  for (int j = 0; j < M; j += 8) {
    for (int ib = 0; ib < B; ib += 8) {
      asm volatile(
        "vlds $39, 0(%[BIAS])\n\t"  "ldi $36, 0(%[IN])\n\t"
        "vlds $32, 0(%[IN])\n\t"    "ldi $37, 0(%[WT])\n\t"
        "vlds $33, 0(%[WT])\n\t"    "ldi $38, 512($31)\n\t"
        // "vsubs $31, $39, $39\n\t"
        "vlenma $39, %[ONE], $31, $40\n\t"
        "vlenma $39, %[ONE], $31, $41\n\t"
        "vlenma $39, %[ONE], $31, $42\n\t"
        "vlenma $39, %[ONE], $31, $43\n\t"
        "vlenma $39, %[ONE], $31, $44\n\t"
        "vlenma $39, %[ONE], $31, $45\n\t"
        "vlenma $39, %[ONE], $31, $46\n\t"
        "vlenma $39, %[ONE], $31, $47\n\t"
        /*"vlds $40, 0*"INSTRIDE"(%[OUT])\n\t"*/ 
        /*"vlds $41, 1*"INSTRIDE"(%[OUT])\n\t"*/ 
        /*"vlds $42, 2*"INSTRIDE"(%[OUT])\n\t"*/ 
        /*"vlds $43, 3*"INSTRIDE"(%[OUT])\n\t"*/
        /*"vlds $44, 4*"INSTRIDE"(%[OUT])\n\t"*/
        /*"vlds $45, 5*"INSTRIDE"(%[OUT])\n\t"*/
        /*"vlds $46, 6*"INSTRIDE"(%[OUT])\n\t"*/
        /*"vlds $47, 7*"INSTRIDE"(%[OUT])\n\t"*/


        "1:\n\t"

        MM_ITER("$32", "$33", "$34", "$35",  " 1*"INSTRIDE"($36)", " 1*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  " 2*"INSTRIDE"($36)", " 2*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  " 3*"INSTRIDE"($36)", " 3*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  " 4*"INSTRIDE"($36)", " 4*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  " 5*"INSTRIDE"($36)", " 5*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  " 6*"INSTRIDE"($36)", " 6*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  " 7*"INSTRIDE"($36)", " 7*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  " 8*"INSTRIDE"($36)", " 8*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  " 9*"INSTRIDE"($36)", " 9*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "10*"INSTRIDE"($36)", "10*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "11*"INSTRIDE"($36)", "11*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "12*"INSTRIDE"($36)", "12*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "13*"INSTRIDE"($36)", "13*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "14*"INSTRIDE"($36)", "14*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "15*"INSTRIDE"($36)", "15*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "16*"INSTRIDE"($36)", "16*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "17*"INSTRIDE"($36)", "17*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "18*"INSTRIDE"($36)", "18*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "19*"INSTRIDE"($36)", "19*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "20*"INSTRIDE"($36)", "20*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "21*"INSTRIDE"($36)", "21*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "22*"INSTRIDE"($36)", "22*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "23*"INSTRIDE"($36)", "23*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "24*"INSTRIDE"($36)", "24*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "25*"INSTRIDE"($36)", "25*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "26*"INSTRIDE"($36)", "26*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "27*"INSTRIDE"($36)", "27*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "28*"INSTRIDE"($36)", "28*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "29*"INSTRIDE"($36)", "29*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "30*"INSTRIDE"($36)", "30*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "31*"INSTRIDE"($36)", "31*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "32*"INSTRIDE"($36)", "32*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "33*"INSTRIDE"($36)", "33*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "34*"INSTRIDE"($36)", "34*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "35*"INSTRIDE"($36)", "35*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "36*"INSTRIDE"($36)", "36*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "37*"INSTRIDE"($36)", "37*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "38*"INSTRIDE"($36)", "38*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "39*"INSTRIDE"($36)", "39*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "40*"INSTRIDE"($36)", "40*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "41*"INSTRIDE"($36)", "41*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "42*"INSTRIDE"($36)", "42*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "43*"INSTRIDE"($36)", "43*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "44*"INSTRIDE"($36)", "44*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "45*"INSTRIDE"($36)", "45*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "46*"INSTRIDE"($36)", "46*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "47*"INSTRIDE"($36)", "47*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "48*"INSTRIDE"($36)", "48*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "49*"INSTRIDE"($36)", "49*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "50*"INSTRIDE"($36)", "50*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "51*"INSTRIDE"($36)", "51*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "52*"INSTRIDE"($36)", "52*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "53*"INSTRIDE"($36)", "53*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "54*"INSTRIDE"($36)", "54*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "55*"INSTRIDE"($36)", "55*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "56*"INSTRIDE"($36)", "56*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "57*"INSTRIDE"($36)", "57*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "58*"INSTRIDE"($36)", "58*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "59*"INSTRIDE"($36)", "59*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "60*"INSTRIDE"($36)", "60*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "61*"INSTRIDE"($36)", "61*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "62*"INSTRIDE"($36)", "62*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$32", "$33", "$34", "$35",  "63*"INSTRIDE"($36)", "63*"WTSTRIDE"($37)", DO) "\n\t"
        MM_ITER("$34", "$35", "$32", "$33",  "64*"INSTRIDE"($36)", "64*"WTSTRIDE"($37)", DO) "\n\t"
        "ldi $36, 64*"INSTRIDE"($36)\n\t"
        "ldi $37, 64*"WTSTRIDE"($37)\n\t"
        "ldi $38, -64($38)\n\t"
        "bne $38, 1b\n\t"

        "vsts $40, 0*"INSTRIDE"(%[OUT])\n\t"
        "vsts $41, 1*"INSTRIDE"(%[OUT])\n\t"
        "vsts $42, 2*"INSTRIDE"(%[OUT])\n\t"
        "vsts $43, 3*"INSTRIDE"(%[OUT])\n\t"
        "vsts $44, 4*"INSTRIDE"(%[OUT])\n\t"
        "vsts $45, 5*"INSTRIDE"(%[OUT])\n\t"
        "vsts $46, 6*"INSTRIDE"(%[OUT])\n\t"
        "vsts $47, 7*"INSTRIDE"(%[OUT])\n\t"
        :: [OUT]"r"(&lout[j][ib]), [IN]"r"(&lin[0][ib]), [WT]"r"(&lnw_tran[0][j]), [BIAS]"r"(bias), [ONE]"r"(simd_vcpyfs(1.0))
        : "$32", "$33", "$34", "$35", "$36", "$37", "$38", "$39",
          "$40", "$41", "$42", "$43", "$44", "$45", "$46", "$47", "memory");

    }
  }

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

__ldm_kernel(linear) __attribute__((aligned(32))) float lout[2][M][B], lin[2][N][B], lnw_tran[N][M], lbias[GM];
__ldm_kernel(linear) volatile long dma_rply = 0, rma_rply = 0, edma_rply, pdma_rply;
static const int x = sizeof(lout) + sizeof(lin) + sizeof(lnw_tran);
void linear(struct mm_par_t *par){
  lwpf_enter(MM);
  lwpf_start(ALL);

  dma_rply = 0;
  dma_iget(par->gbias, lbias, GM*sizeof(float), &dma_rply);
  unsafe_wait_value(&dma_rply, 1);
  for (int j = _MYID * M; j < GM; j += 64*M) {
    dma_rply = 0;
    for (int i = 0; i < GN; i += N) {
      lwpf_start(GET2);
      dma_rply = 0;
      dma_iget(&(par->gw[j/M][i][0]), lnw_tran, N*M*sizeof(float), &dma_rply);
      for (int jj = 0; jj < GM; jj ++)
        lbias[jj] = -lbias[jj];
      unsafe_wait_value(&dma_rply, 1);
      lwpf_stop(GET2);

      edma_rply = 0;
      rma_rply = 0;
      pdma_rply = 0;

      if (_MYID == 0) {
        dma_iget(&(par->gin[0][i][0]), lin[0], N*B*sizeof(float), &edma_rply);
      }
      if (_MYID != 0) {
        asm volatile("sync %0\n\tsynr %0\n\t" :: "r"(0xff) : "memory");
      }
      for (int kk = 0; kk < par->nB; kk += 1) {
        lwpf_start(GET1);
        if (_MYID == 0) {
          unsafe_wait_value(&edma_rply, kk + 1);
          if (kk < par->nB - 1)
            dma_iget(&(par->gin[kk+1][i][0]), lin[kk+1&1], N*B*sizeof(float), &edma_rply);
          // unsafe_wait_value(&rma_rply, kk);
          asm volatile("sync %0\n\tsynr %0\n\t" :: "r"(0xff) : "memory");
          rma(0, lin[kk&1], lin[kk&1], N*B*sizeof(float), &rma_rply, &rma_rply, RMA_MODE_COL_BCAST, RMA_OPCODE_PUT, 0xff);
        } else {
          unsafe_wait_value(&rma_rply, kk+1);
          asm volatile("sync %0\n\tsynr %0\n\t" :: "r"(0xff) : "memory");
        }
        lwpf_stop(GET1);
        lwpf_start(SIMD);
        plain_matmul(lout[kk&1], lin[kk&1], lnw_tran, lbias + j);
        lwpf_stop(SIMD);
        lwpf_start(PUT);
        unsafe_wait_value(&pdma_rply, kk);
        // dma_rply = 0;
        dma_iput(&(par->gout[kk][j][0]), lout[kk&1], M*B*sizeof(float), &pdma_rply);
        // unsafe_wait_value(&dma_rply, 1);
        lwpf_stop(PUT);
      }
      if (_MYID == 0) {
        asm volatile("sync %0\n\tsynr %0\n\t" :: "r"(0xff) : "memory");
      }
      unsafe_wait_value(&pdma_rply, par->nB);
    }
  }
  lwpf_stop(ALL);
  // }
  lwpf_exit(MM);
}
#pragma omp declare target(linear)
#endif

#ifdef __sw_host__
#ifdef MATMUL_TEST
__attribute__((aligned(256))) float goutpp[GB/B][GM][B], ginpp[GB/B][GN][B], gwpp[GM/M][GN][M];
float gout[GB][GM], gin[GB][GN], gw[GM][GN], gout_ref[GB][GN], gbias[GM];
void init(){
  for (int ib = 0; ib < GB; ib ++) {
    for (int i = 0; i < GN; i ++){
      gin[ib][i] = rand() * 1.0/ INT_MAX;
    }
  }
  for (int i = 0; i < GN; i ++){
    for (int j = 0; j < GM; j ++) {
      gw[i][j] = rand() * 1.0/ INT_MAX;
    }
  }
  for (int j = 0; j < GM; j ++)
    gbias[j] = rand() * 1.0/ INT_MAX;
  for (int ib = 0; ib < GB; ib++) {
    for (int j = 0; j < GM; j++) {
      gout_ref[ib][j] = gbias[j];
      for (int i = 0; i < GN; i++) {
        gout_ref[ib][j] += gin[ib][i] * gw[j][i];
      }
    }
  }
}
void convert(){
  for (int i = 0; i < GN; i ++) {
    for (int j = 0; j < GM; j ++) {
      gwpp[j/M][i][j%M] = -gw[j][i];
    }
  }
  for (int ib = 0; ib < GB; ib ++) {
    for (int i = 0; i < GN; i ++) {
      goutpp[ib / B][i][ib % B] = gout[ib][i];
    }
    for (int j = 0; j < GM; j ++) {
      ginpp[ib/B][j][ib%B] = gin[ib][j];
    }
  }
}
void convertb(){
  for (int ib = 0; ib < GB; ib ++) {
    for (int i = 0; i < GN; i ++) {
      gout[ib][i] = goutpp[ib/B][i][ib%B];
    }
  }
}

#include <athread.h>
#define LWPF_UNITS U(MM)
#include "lwpf3/lwpf.h"
extern void slave_linear(void *);
int main(){
  init();
  convert();
  // linear(gin, gw, gout_ref, )
  athread_init();
  lwpf_init(NULL);
  struct mm_par_t par = {goutpp, ginpp, gwpp, gbias, (GB+B-1) / B};
  __real_athread_spawn(slave_linear, &par, 1);
  athread_join();
  convertb();
  float abserr = 0;
  for (int i = 0; i < GB; i ++)
    for (int j = 0; j < GM; j ++) {
      if (fabs(gout[i][j] - gout_ref[i][j]) > abserr) {
        abserr = fabs(gout[i][j] - gout_ref[i][j]);
        printf("%d %d %f %f %f\n", i, j, gout_ref[i][j], gout[i][j], abserr);
      }
    }
  lwpf_report_summary(stdout);
  lwpf_init(NULL);
  athread_enter64_arg();
  __real_athread_spawn64_arg(slave_linear, &par);
  athread_join64_arg();
  athread_leave64_arg();
  lwpf_report_summary(stdout);
}
#endif
#endif