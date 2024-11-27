#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <assert.h>
extern void *__real_malloc(size_t size);
extern void __real_free(void*);
extern void *__real_calloc(size_t n, size_t size);
extern void *__real_realloc(void *ptr, size_t size);
const size_t POOL_BASE_GROWTH = (64*1024UL*1024UL);
struct segment {
  void *base;
  size_t size;
  struct segment *next;
};
struct mempool{
  struct segment *first;
  size_t block_size;
  void **free_ptrs;
  int nfree;
  int cap;
  int nall;
};
static void init_pool(struct mempool *pool, size_t size){
  pool->first = NULL;
  pool->block_size = size;
  pool->cap = POOL_BASE_GROWTH/size < 16384 ? 16384 : POOL_BASE_GROWTH/size;
  pool->free_ptrs = __real_malloc(pool->cap * sizeof(void*));
  pool->nfree = 0;
  pool->nall = 0;
}
static int busy = 0;
static int ngrow = 0;
static void grow_free_ptrs(struct mempool *pool){
  void **new_free_ptrs = __real_malloc(pool->cap * 2 * sizeof(void*));
  for (int i = 0; i < pool->nfree; i ++)
    new_free_ptrs[i] = pool->free_ptrs[i];
  __real_free(pool->free_ptrs);
  pool->free_ptrs = new_free_ptrs;
  pool->cap *= 2;
}
static void grow_pool(struct mempool *pool){
  ngrow ++;
  // printf("%ld\n", POOL_BASE_GROWTH + sizeof(struct segment) + 256);
  struct segment *new_segment = __real_malloc(POOL_BASE_GROWTH + sizeof(struct segment) + 256);
  // printf("%p\n", new_segment);
  assert(new_segment); 
  new_segment->base = (void*)(((intptr_t)new_segment + sizeof(struct segment) + 255L) / 256L * 256L);
  new_segment->size = POOL_BASE_GROWTH;
  new_segment->next = pool->first;
  pool->first = new_segment;
  size_t nnew = new_segment->size / pool->block_size;
  // busy = 1;
  // printf("%d %p %d %ld %d\n", pool->cap, new_segment, pool->nfree, pool->block_size, pool->nall);
  // busy = 0;

  while (pool->nall + nnew > pool->cap) {
    grow_free_ptrs(pool);
  }

  for (size_t i = 0; i < nnew; i ++) {
    void *blk = (void*)((intptr_t)new_segment->base + (pool->block_size * i));
    pool->free_ptrs[pool->nfree++] = blk;
  }
  pool->nall += nnew;
}
static long nsearch = 0, niter = 0;
static int ptr_in_pool(void *ptr, struct mempool *pool){
  for (struct segment *seg = pool->first; seg; seg = seg->next) {
    niter ++;
    if (ptr >= seg->base && ptr < seg->base + seg->size) 
      return 1;
  }
  return 0;
}

static const size_t pool_sizes[] = {
  256,
  1024-256,
  2 * 1024,
  4 * 1024,
  8 * 1024,
  16 * 1024,
  32 * 1024-256,
  64 * 1024-256,
  128 * 1024-256,
  256 * 1024-256,
  512 * 1024-768,
  1024 * 1024-256,
  2 * 1024 * 1024-256,
  4 * 1024 * 1024-256,
  8 * 1024 * 1024-256,
  16 * 1024 * 1024-256,
  32 * 1024 * 1024-256
};
#define NPOOL (sizeof(pool_sizes)/sizeof(pool_sizes[0]))
static struct mempool pools[sizeof(pool_sizes)/sizeof(pool_sizes[0])]; //64, 256, 1K, 2K, 4K, 8K, 16K, 32K, 64K, 128K, 256K, 512K, 1M, 2M, 4M, 8M, 16M
#define GUARD_POOL (pools + NPOOL)
static struct mempool *find_ptr_pool(void *ptr){
  nsearch ++;
  for (int i = 0; i < NPOOL; i ++) {
    if (ptr_in_pool(ptr, pools+i)) {
      return pools + i;
    }
  }
  return NULL;
}
void mempools_report_(){
  size_t tot_size = 0;
  for (int i = 0; i < NPOOL; i ++){
    tot_size += pools[i].block_size * pools[i].nall;
    printf("number of block size=%ld is %d\n", pools[i].block_size, pools[i].nall - pools[i].nfree);
  }
  printf("total pool size: %ld, num searches=%ld, actually accessed segment=%ld\n", tot_size, nsearch, niter);

}

static struct mempool *find_size_pool(size_t sz){
  struct mempool *pool = pools;
  while (pool < GUARD_POOL && pool->block_size < sz) pool ++;
  if (pool == GUARD_POOL) return NULL;
  return pool;
}
static int pool_initialized = 0;
int pool_enable = 0;
void mempools_init_(){
  if (pool_initialized) return;
  for (int i = 0; i < NPOOL; i ++) {
    init_pool(pools + i, pool_sizes[i]);
  }
  pool_initialized = 1;
}
void mempools_enable_(){
  pool_enable = 1;
}
void mempools_disable_(){
  pool_enable = 0;
}
static void *get_block(struct mempool *pool) {
  if (pool->nfree == 0) {
    grow_pool(pool);
  }
  return pool->free_ptrs[--pool->nfree];
}
void *pool_malloc(size_t sz){
  // if (busy) 
  // return __real_malloc(sz);
  if (!pool_enable || busy) return __real_malloc(sz);
  struct mempool *pool = find_size_pool(sz);
  if (pool) {
    void *ret = get_block(pool);
    // if (((intptr_t)ret & 0xfff000000000L) != 0x500000000000L || !ptr_in_pool(ret, pool)){
    //   busy = 1;
    //   printf("strange ptr: %p %d\n", ret, pool - pools);
    //   busy = 0;
    // }
    // memset(ret, 0, sz);
    return ret;
  } else {
    void *ret = __real_malloc(sz);
    return ret;
  }
}
void pool_free(void *ptr){
  // __real_free(ptr);
  // return;
  if (!pool_initialized || busy) {
    __real_free(ptr);
    return;
  }
  struct mempool *pool = find_ptr_pool(ptr);
  if (pool) {
    if (((intptr_t)ptr & 0xfff000000000L) != 0x500000000000L){
      busy = 1;
      printf("strange ptr: %p %d\n", ptr, pool - pools);
      busy = 0;
    }
    pool->free_ptrs[pool->nfree++] = ptr;
  }
  else
    __real_free(ptr);
}
void *pool_realloc(void *ptr, size_t sz){
  if (!pool_enable || busy) return __real_realloc(ptr, sz);
  struct mempool *pool = find_ptr_pool(ptr);
  if (pool) {
    if (pool->block_size >= sz) {
      return ptr;
    } else {
      void *ret = pool_malloc(sz);
      memcpy(ret, ptr, pool->block_size);
      // memset(ret + pool->block_size, 0, sz - pool->block_size);
      if (((intptr_t)ptr & 0xfff000000000L) != 0x500000000000L || !ptr_in_pool(ptr, pool)){
        busy = 1;
        printf("strange ptr: %p %d\n", ptr, pool - pools);
        busy = 0;
      }
      pool->free_ptrs[pool->nfree++] = ptr;
      return ret;
    }
  } else {
    return __real_realloc(ptr, sz);
  }
}
void *pool_calloc(size_t n, size_t sz) {
  if (!pool_enable || busy) return __real_calloc(n, sz);
  void *ret = pool_malloc(sz*n);
  // busy = 1;
  memset(ret, 0, sz*n);
  // busy = 0;
  // busy = 1;
  // printf("%p %d\n", ret, sz*n);
  // busy = 0;
  return ret;
}
void *__wrap_malloc(size_t) __attribute__((alias("pool_malloc")));
void *__wrap_calloc(size_t, size_t) __attribute__((alias("pool_calloc")));
void *__wrap_realloc(void *, size_t) __attribute__((alias("pool_realloc")));
void __wrap_free(void *) __attribute__((alias("pool_free")));
// #include <mpi.h>
// int main(int argc, char **argv){
//   puts("!!!");
//   // sleep(10);
//   MPI_Init(&argc, &argv);
// }
// FILE *mlog = NULL;
// int trace = 1;
// void trace_malloc_start_(int *rank){
//   char fname[12] = "mlog.";
//   char *p = fname + 5;
//   *(p++) = *rank / 10000 + '0';
//   *(p++) = *rank % 10000 / 1000 + '0';
//   *(p++) = *rank % 1000 / 100 + '0';
//   *(p++) = *rank % 100 / 10 + '0';
//   *(p++) = *rank % 10 + '0';
//   *(p++) = '\0';
//   mlog = fopen(fname, "w"); //open(fname, O_RDWR|O_CREAT|O_TRUNC, 0644);
// }
// void *__wrap_malloc(size_t size){
//   void *ret = __real_malloc(size);
//   if (mlog && trace) {
//     trace = 0;
//     fprintf(mlog, "malloc %p %ld\n", ret, size);
//     trace = 1;
//   }
//   return ret;
// }
// void __wrap_free(void *ptr) {
//   __real_free(ptr);
//   if (mlog && trace) {
//     trace = 0;
//     fprintf(mlog, "free %p\n", ptr);
//     trace = 1;
//   }
// }
// void *__wrap_calloc(size_t n, size_t size) {
//   void *ret = __real_calloc(n, size);
//   if (mlog && trace) {
//     trace = 0;
//     fprintf(mlog, "calloc %p %ld\n", ret, n*size);
//     trace = 1;
//   }
//   return ret;
// }
// void *__wrap_realloc(void *ptr, size_t size) {
//   void *ret = __real_realloc(ptr, size);
//   if (mlog && trace) {
//     trace = 0;
//     fprintf(mlog, "realloc %p %p %ld\n", ptr, ret, size);
//     trace = 1;
//   }
//   return ret;
// }