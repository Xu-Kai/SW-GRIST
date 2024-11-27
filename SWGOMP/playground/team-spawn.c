extern int swgomp_get_cpe_num();
#include <omp.h>
int main(){
  #pragma omp target teams num_teams(3)
  {
    cal_locked_printf("team: %d %d\n", swgomp_get_cpe_num(), omp_get_team_num());
    if (omp_get_team_num() == 2) {
      #pragma omp parallel num_threads(4)
      {
        cal_locked_printf("thread: %d %d %d\n", swgomp_get_cpe_num(), omp_get_team_num(), omp_get_thread_num());
      }
    }
    if (omp_get_team_num() == 1) {
      omp_set_nested(1);
      #pragma omp parallel num_threads(4)
      {
        cal_locked_printf("thread: %d %d %d\n", swgomp_get_cpe_num(), omp_get_team_num(), omp_get_thread_num());
      }
    }
  }
}