#ifdef __sw_slave__
#define RFPCR0(x) {asm volatile("rfpcr $63\n\tvextf $63,0,%0":"=r"(x)::);}
#define RFPCR1(x) {asm volatile("rfpcr $63\n\tvextf $63,1,%0":"=r"(x)::);}
#define RFPCR(x,y) {asm volatile("rfpcr $63\n\tvextf $63,0,%0\n\tvextf $63,1,%1":"=r"(x),"=r"(y)::);}
#define WFPCR(x,y) {asm volatile("vinsf %0,$63,0,$63\n\tvinsf %1,$63,1,$63\n\twfpcr $63"::"r"(x),"r"(y):"memory");}
int is_ftz() {
  long v = 0, vh;
  RFPCR0(v);
  RFPCR1(vh);
  return (v & (1L << 48)) == 0;
}
void enable_ftz(){
  long v = 0, vh;
  RFPCR0(v);
  RFPCR1(vh);
  v &= ~(1L << 48);
  WFPCR(v, vh);
}
void disable_ftz(){
  long v = 0, vh;
  RFPCR0(v);
  RFPCR1(vh);
  v |= (1L << 48);
  WFPCR(v, vh);
}
#endif
#ifdef __sw_host__
#include <stddef.h>
#include "athread_compat.h"
extern void slave_enable_ftz();
extern void slave_disable_ftz();
void enable_ftz_all_cpe(){
  jobserver_athread_spawn(slave_enable_ftz, NULL, 1);
  jobserver_athread_join();
}
void disable_ftz_all_cpe(){
  jobserver_athread_spawn(slave_enable_ftz, NULL, 1);
  jobserver_athread_join();
}
void enable_ftz_all_cpe_() __attribute__((alias("enable_ftz_all_cpe")));
void disable_ftz_all_cpe_() __attribute__((alias("disable_ftz_all_cpe")));
#endif