int add(int n, int *a, int *b, int *c){
  #pragma omp target parallel for schedule(static,64)
  for (int i = 0; i < n; i ++) {
    a[i] = b[i] + c[i];
  }
}