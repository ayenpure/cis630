//pagerankmpi
#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  int num_tasks, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cout << "Process " << rank << " of " << num_tasks << endl;
  MPI_Finalize();
  return 0;
}
