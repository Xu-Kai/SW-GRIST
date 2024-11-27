#include <stdint.h>
#include <stdio.h>
static uint64_t calc_bitmap(int rootpe, int clev, int nthr) {
  uint64_t bm = 0;
  struct {
    int root;
    int clev;
    int nthr;
  } q[64];
  q[0] = (typeof(q[0])){rootpe, clev, nthr};
  int s = 0, t = 1;
  while (s < t) {
    int root = q[s].root;
    int clev = q[s].clev;
    int nthr = q[s].nthr;
    int nl = (nthr + 1) >> 1;
    int nr = nthr - nl;
    if (nr > 0) {
      q[t++] = (typeof(q[0])){root + (1 << clev), clev - 1, nr};
    }
    if (nl > 1) {
      q[s] = (typeof(q[0])){root, clev - 1, nl};
    } else {
      bm |= (uint64_t)1 << root;
      s++;
    }
  }
  return bm;
}
int main(){
  uint64_t bm = calc_bitmap(0, 5, 6);
  for (int i = 0; i < 8; i ++) {
    for (int j = 0; j < 8; j ++) {
      int id = i * 8 + j;
      printf("%d", (bm & 1UL << id) != 0);
    }
    puts("");
  }
}