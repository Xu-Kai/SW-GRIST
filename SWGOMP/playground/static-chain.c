// #pragma omp declare target
int test(){
  const char *s = "abc";
  void exec(int n) {
    for (int i = 0; i < 6; i ++)
      puts(s);
  }
  exec(6);
}
// #pragma omp end declare target
#pragma swgomp entries(test, exec)
int x;