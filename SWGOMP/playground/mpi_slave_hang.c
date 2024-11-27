#include <mpi.h>

int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  swgomp_set_slave_wait_signal(44);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #pragma omp target parallel
  {
    if (rank == 11) while (1);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}