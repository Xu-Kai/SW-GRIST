#ifdef __sw_slave__
#include <stdint.h>
#include "swdefs.h"
const intptr_t PTR_UNCACHED_BIT = 0x200000000000L;
void *uncached_ptr(void *ptr) {
  return (void*)((intptr_t)ptr | PTR_UNCACHED_BIT);
}
#define MAX(a, b) (((a)>(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))
#define DEF_CAS(NAME, TYPE, SFI, SFS, CMP, SWP)                                \
  static inline long NAME(TYPE *ptr, TYPE cmp, TYPE swp) {                     \
    int512 tmp;                                                                \
    asm volatile("vins" #SFI " %2, %0, " #CMP ", %0\n\t"                       \
                 "vins" #SFI " %3, %0, " #SWP ", %0\n\t"                       \
                 "cas"  #SFS " %0, 0(%4)\n\t"                                  \
                 : "=&r"(tmp), "+m"(*ptr)                                      \
                 : "r"(cmp), "r"(swp), "r"(ptr));                              \
    return *(long *)&tmp;                                                      \
  }
DEF_CAS(casw , int, w, w, 2, 4)
DEF_CAS(caswu, unsigned int, w, w, 2, 4)
union fw {
  float f;
  int w;
};
static inline long casf(float *ptr, float cmp, float swp) {
  int512 tmp;
  union fw ucmp = {.f = cmp}, uswp = {.f = swp};
  asm volatile("vinsw %2, %0, 2, %0\n\t"
               "vinsw %3, %0, 4, %0\n\t"
               "casw  %0, 0(%4)\n\t"
               : "=&r"(tmp), "+m"(*ptr)
               : "r"(ucmp.w), "r"(uswp.w), "r"(ptr));
  return *(long *)&tmp;
}
DEF_CAS(casl , long, f, l, 1, 2)
DEF_CAS(caslu, unsigned long, f, l, 1, 2)
DEF_CAS(casd , double, f, l, 1, 2)

#define DEF_REDUCE(OP, TYPE, MODE, CAS, EXPR)                                  \
  void swgomp_reduce_##OP##_##MODE(TYPE *orig_ptr, TYPE val) {                 \
    TYPE *orig_ptr_uc = uncached_ptr(orig_ptr);                                \
    long success;                                                              \
    do {                                                                       \
      TYPE orig_val = *orig_ptr_uc;                                            \
      TYPE new_val = EXPR;                                                     \
      success = CAS(orig_ptr_uc, orig_val, new_val);                           \
    } while (!success);                                                        \
  }

DEF_REDUCE(plus       , int, SI, casw, orig_val + val)
DEF_REDUCE(mult       , int, SI, casw, orig_val * val)
DEF_REDUCE(max        , int, SI, casw, MAX(orig_val, val))
DEF_REDUCE(min        , int, SI, casw, MIN(orig_val, val))
DEF_REDUCE(truth_andif, int, SI, casw, orig_val && val)
DEF_REDUCE(truth_orif , int, SI, casw, orig_val || val)
DEF_REDUCE(bit_and    , int, SI, casw, orig_val & val)
DEF_REDUCE(bit_or     , int, SI, casw, orig_val | val)
DEF_REDUCE(bit_xor    , int, SI, casw, orig_val ^ val)

DEF_REDUCE(plus       , unsigned int, SI_unsigned, caswu, orig_val + val)
DEF_REDUCE(mult       , unsigned int, SI_unsigned, caswu, orig_val * val)
DEF_REDUCE(max        , unsigned int, SI_unsigned, caswu, MAX(orig_val, val))
DEF_REDUCE(min        , unsigned int, SI_unsigned, caswu, MIN(orig_val, val))
DEF_REDUCE(truth_andif, unsigned int, SI_unsigned, caswu, orig_val && val)
DEF_REDUCE(truth_orif , unsigned int, SI_unsigned, caswu, orig_val || val)
DEF_REDUCE(bit_and    , unsigned int, SI_unsigned, caswu, orig_val & val)
DEF_REDUCE(bit_or     , unsigned int, SI_unsigned, caswu, orig_val | val)
DEF_REDUCE(bit_xor    , unsigned int, SI_unsigned, caswu, orig_val ^ val)

DEF_REDUCE(plus       , long, DI, casl, orig_val + val)
DEF_REDUCE(mult       , long, DI, casl, orig_val * val)
DEF_REDUCE(max        , long, DI, casl, MAX(orig_val, val))
DEF_REDUCE(min        , long, DI, casl, MIN(orig_val, val))
DEF_REDUCE(truth_andif, long, DI, casl, orig_val && val)
DEF_REDUCE(truth_orif , long, DI, casl, orig_val || val)
DEF_REDUCE(bit_and    , long, DI, casl, orig_val & val)
DEF_REDUCE(bit_or     , long, DI, casl, orig_val | val)
DEF_REDUCE(bit_xor    , long, DI, casl, orig_val ^ val)

DEF_REDUCE(plus       , unsigned long, DI_unsigned, caslu, orig_val + val)
DEF_REDUCE(mult       , unsigned long, DI_unsigned, caslu, orig_val * val)
DEF_REDUCE(max        , unsigned long, DI_unsigned, caslu, MAX(orig_val, val))
DEF_REDUCE(min        , unsigned long, DI_unsigned, caslu, MIN(orig_val, val))
DEF_REDUCE(truth_andif, unsigned long, DI_unsigned, caslu, orig_val && val)
DEF_REDUCE(truth_orif , unsigned long, DI_unsigned, caslu, orig_val || val)
DEF_REDUCE(bit_and    , unsigned long, DI_unsigned, caslu, orig_val & val)
DEF_REDUCE(bit_or     , unsigned long, DI_unsigned, caslu, orig_val | val)
DEF_REDUCE(bit_xor    , unsigned long, DI_unsigned, caslu, orig_val ^ val)

DEF_REDUCE(plus       , float, SF, casf, orig_val + val)
DEF_REDUCE(mult       , float, SF, casf, orig_val * val)
DEF_REDUCE(max        , float, SF, casf, MAX(orig_val, val))
DEF_REDUCE(min        , float, SF, casf, MIN(orig_val, val))
DEF_REDUCE(truth_andif, float, SF, casf, orig_val && val)
DEF_REDUCE(truth_orif , float, SF, casf, orig_val || val)

DEF_REDUCE(plus       , double, DF, casd, orig_val + val)
DEF_REDUCE(mult       , double, DF, casd, orig_val * val)
DEF_REDUCE(max        , double, DF, casd, MAX(orig_val, val))
DEF_REDUCE(min        , double, DF, casd, MIN(orig_val, val))
DEF_REDUCE(truth_andif, double, DF, casd, orig_val && val)
DEF_REDUCE(truth_orif , double, DF, casd, orig_val || val)
// void swgomp_reduce_max_SI(int *orig_ptr, int val) {
//   int *orig_ptr_uc = uncached_ptr(orig_ptr);
//   long success = 0;
//   while (!success) {
//     int orig_val = *orig_ptr_uc;
//     int new_val = MAX(orig_val, val);
//     success = casw(orig_ptr_uc, orig_val, new_val);
//   }
// }

#endif