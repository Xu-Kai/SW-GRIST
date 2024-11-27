/* This is suggestted for SW5 since they do not have __sw_{host,slave}__ defs*/
#include "lwpf_where.h"
#ifdef __sw_host__
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <athread.h>
#define LWPF_UNITS U(TEST)
#include "lwpf.h"
extern void slave_lwpf_test();
int main(){
  athread_init();
  lwpf_init(NULL);
  athread_spawn(lwpf_test, 0);
  athread_join();
  lwpf_report_summary(stdout);
  //lwpf_report_detail(stdout);

  evt_conf_t conf;
#if defined(SW5)
  conf.pc_mask = 0xf;
  conf.evt[0] = PC0_CNT_INST;
  conf.evt[1] = PC1_CYCLE;
  conf.evt[2] = PC2_CNT_GST;
  conf.evt[3] = PC3_CYCLE;
#elif defined(SW7)
  conf.pc_mask = 0xff;
  conf.evt[0] = PC0_CYCLE;
  conf.evt[1] = PC1_INST;
  conf.evt[2] = PC2_L1IC_ACCESS;
  conf.evt[3] = PC3_GST;
  conf.evt[4] = PC4_L1IC_MISSTIME;
  conf.evt[5] = PC5_INST_L0IC_READ;
  conf.evt[6] = PC6_CYC_LAUNHCNONE_BUFFER;
  conf.evt[7] = PC7_INST_L0IC_READ;
#elif defined(SW9)
  conf.pc_mask = 0xff;
  conf.evt[0] = PC0_CYCLE;
  conf.evt[1] = PC1_INST;
  conf.evt[2] = PC2_L1IC_ACCESS;
  conf.evt[3] = PC3_GST;
  conf.evt[4] = PC4_L1IC_MISSTIME;
  conf.evt[5] = PC5_INST_L0IC_READ;
  conf.evt[6] = PC6_CYC_LAUNHCNONE_BUFFER;
  conf.evt[7] = PC7_INST_L0IC_READ;
#endif
  lwpf_init(&conf);
  athread_spawn(lwpf_test, 0);
  athread_join();
  lwpf_report_summary(stdout);
  //lwpf_report_detail(stdout);

}
#endif
#ifdef __sw_slave__
#include <simd.h>
#include <slave.h>
#if defined(SW7) || defined(SW9)
#include <crts.h>
#endif
#include <stdio.h>
#if defined(SW7) || defined(SW9)
#define EVT_PC0 PC0_CYCLE                
#define EVT_PC1 PC1_INST                 
#define EVT_PC2 PC2_L1IC_ACCESS          
#define EVT_PC3 PC3_LDM_PIPE             
#define EVT_PC4 PC4_L1IC_MISSTIME        
#define EVT_PC5 PC5_INST_L0IC_READ       
#define EVT_PC6 PC6_CYC_LAUNHCNONE_BUFFER
#define EVT_PC7 PC7_INST_L0IC_READ
#elif defined(SW5)
#define EVT_PC0 PC0_CNT_INST
#define EVT_PC1 PC1_CYCLE
#define EVT_PC2 PC2_CNT_GLD
#define EVT_PC3 PC3_CYCLE
#endif

#define LWPF_KERNELS K(A) K(B) K(C)
#define LWPF_UNIT U(TEST)
#include "lwpf.h"
void lwpf_test(){
  lwpf_enter(TEST);
  lwpf_start(C);
  lwpf_stop(C);
  lwpf_start(A);
  lwpf_start(B);
  lwpf_stop(B);
  lwpf_stop(A);
  lwpf_exit(TEST);
}
#endif
