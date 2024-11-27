#include <string.h>
#include "swdefs.h"
__always_inline int in_ldm(void *ptr){
  return (ptr > (void*)0 && ptr < (void*)0x40000);
}
__always_inline void dma_get(void *mem, void *ldm, size_t count) {
  //long dma_desc[4];
  int512 dma_desc;
  asm volatile(
      "ble %[SIZE], 2f\n\t"
      "memb\n\t"
      "stl $31, 0(%[RPL])\n\t"
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "ldi %[DESC], 1($31)\n\t"
      "sll %[DESC], 56, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "1:\n\t"
      "ldw %[DESC], 0(%[RPL])\n\t"
      "beq %[DESC], 1b\n\t"
      "2:\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(((long*)&dma_desc) + 2), [ SIZE ] "r"(count)
      : "memory");
}
__always_inline void dma_put(void *mem, void *ldm, size_t count) {
  int512 dma_desc;
  asm volatile(
      "ble %[SIZE], 2f\n\t"
      "memb\n\t"
      "stl $31, 0(%[RPL])\n\t"
      "vinsw %[RPL], $31, 2, %[DESC]\n\t"
      "vinsw %[SIZE], %[DESC], 0, %[DESC]\n\t"
      "dma %[DESC], %[MEM], %[LDM]\n\t"
      "1:\n\t"
      "ldw %[DESC], 0(%[RPL])\n\t"
      "beq %[DESC], 1b\n\t"
      "2:\n\t"
      "memb\n\t"
      : [ DESC ] "=&r"(dma_desc)
      : [ MEM ] "r"(mem), [ LDM ] "r"(ldm), [ RPL ] "r"(((long*)&dma_desc) + 2), [ SIZE ] "r"(count)
      : "memory");
}

void omnicopy(void *dst, void *src, size_t size) {
#ifndef __sw_slave__
  memcpy(dst, src, size);
#else
  switch (in_ldm(dst) << 1 | in_ldm(src)) {
    case 0: /*fallthrou*//*都不在LDM*/
    case 3: /*都在LDM*/
    memcpy(dst, src, size);
    break;
    case 1: /*src 在LDM*/
    dma_put(dst, src, size);
    break;
    case 2: /*dst 在LDM*/
    dma_get(src, dst, size);
  }
#endif
}

void omnicopy_c_(void *dst, void *src, int *size) {
  // cal_locked_printf("%p %p %ld\n", dst, src, *size);
  omnicopy(dst, src, *size);
}