#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank, num_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  //cout << "Hello World form process " << rank << " of " << num_proc;
  printf("Hello World from process %d of %d\n", rank, num_proc);
}
