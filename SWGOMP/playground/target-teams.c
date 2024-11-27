int add_target(int n, int *a, int *b, int *c){
  #pragma omp target
  for (int i = 0; i < n; i ++) {
    a[i] = b[i] + c[i];
  }
}
int add_teams(int n, int *a, int *b, int *c){
  #pragma omp target teams
  #pragma omp distribute dist_schedule(guided,64)
  for (int i = 0; i < n; i ++) {
    a[i] = b[i] + c[i];
  }
}