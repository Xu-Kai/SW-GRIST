
#include <stdio.h>
#include <omp.h>

extern void *__real_memcpy(void *dest, void *src, size_t size);
extern void omnicopy(void *dst, void *src, size_t size);
inline int pointer_check(void *src){
    return (((long)src & 3)== 0) && (long)src >= 0x500000000000L && (long)src <= 0x600000000000L;
}
inline int alias_check(void *dest, void *src, int size) {
    long ldst = (long)dest;
    long lsrc = (long)src;
    long laddr = ldst < lsrc ? ldst : lsrc;
    long raddr = ldst > lsrc ? ldst : lsrc;
    return laddr + size < raddr;
}
void hybrid_memcpy(void *dest, void *src, size_t size)
{


    //  __real_memcpy(dest, src, size);
    //  return;
    if (size < 64 * 1024 || !pointer_check(src) || !pointer_check(dest) || !alias_check(dest, src, size))
    {
        __real_memcpy(dest, src, size);
    }
    else{
        #pragma omp target parallel for
        for (long i = 0; i < size; i += 2048) {
            int buf[512];
            int r = size - i < 2048 ? size - i : 2048;
            omnicopy(buf, src + i, r);
            omnicopy(dest + i, buf, r);
        }
    }

}
void __wrap_memcpy(void *dest, void *src, size_t size) __attribute__((alias("hybrid_memcpy")));
