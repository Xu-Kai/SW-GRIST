#ifdef __sw_host__
#include <stdio.h>
#define def_lwpf_report(base)                  \
  void base ## _(int *who){                    \
    if (*who == -1)                            \
      base(stdout);                            \
    else{                                      \
      char fn[64];                             \
      FILE *f = fopen(fn,"w");                 \
      base(f);                                 \
    }                                          \
  }
def_lwpf_report(lwpf_report_summary);
def_lwpf_report(lwpf_report_detail);
__attribute__((weak)) void lwpf_init_(evt_conf_t *evtconf){
  //printf("%p\n",evtconf);
  lwpf_init(evtconf);
}

__attribute__((weak)) void lwpf_init_null_() {
  lwpf_init(NULL);
}

#endif

#ifdef __sw_slave__
#define FCAT_INNER(x,y) x ## _ ## y ## _
#define FCAT(x,y) FCAT_INNER(x,y)
#define K(x) FCAT(LWPF_UNIT,x)(){lwpf_start(x);}
#define U(x) void lwpf_stat_ ## x
LWPF_KERNELS
#undef U
#undef K

#define K(x) FCAT(LWPF_UNIT,x)(){lwpf_stop(x);}
#define U(x) void lwpf_stop_ ## x
LWPF_KERNELS
#undef U
#undef K

#define U(x) lwpf_enter_ ## x ## _ () {lwpf_sync_counters_m2c(lwpf_global_counter_ ## x[_MYID],lwpf_kernel_count_ ## x);}
LWPF_UNIT
#undef U
#define U(x) lwpf_exit_ ## x ## _() {lwpf_sync_counters_c2m(lwpf_global_counter_ ## x[_MYID],lwpf_kernel_count_ ## x);}
LWPF_UNIT
#undef U
#undef FCAT
#undef FCAT_INNER
#endif
