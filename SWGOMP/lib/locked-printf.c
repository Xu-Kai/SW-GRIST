#ifdef __sw_slave__
static volatile __uncached long req = 0, cur = 0;
#include <stdio.h>
#include <stdarg.h>
void cal_locked_printf(const char *fmt, ...){
  int ireq, icur;
  asm volatile(
    "faal %[IREQ], 0(%[MREQ])\n\t"
    "1:\n\t"
    "ldl %[ICUR], %[MCUR]\n\t"
    "subl %[ICUR], %[IREQ], %[ICUR]\n\t"
    "bne %[ICUR], 1b\n\t"
  :[IREQ]"=&r"(ireq), [ICUR]"=&r"(icur), "+m"(req), [MCUR]"+m"(cur) : [MREQ]"r"(&req): "memory");
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  asm volatile("faal $31, 0(%1)": "+m"(cur):"r"(&cur):"memory");
}
#endif
