#include <stdlib.h>
#include <stdio.h>
#include <gptl.h>
extern void (*SWGOMP_pre_target_cb)(int dev, void *fn, int b, void *dat, void *sz,
                        void *kinds, int c, size_t d, void **args);
extern void (*SWGOMP_post_target_cb)(int dev, void *fn, int b, void *dat, void *sz,
                        void *kinds, int c, size_t d, void **args);
void pre_target(int dev, void *fn, int b, void *dat, void *sz, void *kinds, int c, size_t d, void **args) {
  char s[64];
  sprintf(s, "func@%p", fn);
  GPTLstart(s);
}
void post_target(int dev, void *fn, int b, void *dat, void *sz, void *kinds, int c, size_t d, void **args) {
  char s[64];
  sprintf(s, "func@%p", fn);
  GPTLstop(s);
}
void gptl_target_enable_(){
  SWGOMP_pre_target_cb = pre_target;
  SWGOMP_post_target_cb = post_target;
}