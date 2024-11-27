#include <omp.h>
void callee(){

}
int main(){
  // #pragma omp target parallel for schedule(dynamic, 16)
  // for (int i = 0; i < 4096; i += 2) {
  //   cal_locked_printf("%d %d %d %d\n", swgomp_get_cpe_num(), omp_get_thread_num(), omp_get_num_threads(), i);
  // }
  // int n;
#pragma omp target
  {
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 4096; i += 2) {
// #pragma omp single
      callee();
      cal_locked_printf("%d %d %d %d\n", swgomp_get_cpe_num(), omp_get_thread_num(), omp_get_num_threads(), i);
    }
  }
}
