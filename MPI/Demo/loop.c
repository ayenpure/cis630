/*
  "Hello World" MPI Test Program
*/
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    char buf[256];
    int my_rank, num_procs;
    /* Initialize the infrastructure necessary for communication */
    MPI_Init(&argc, &argv);
    /* Identify this process */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf("The rank for this process is %d\n", my_rank);
    /* Find out how many total processes are active */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    printf("Out of total %d processes\n", num_procs);

    for(int i = 0; i < 100000; i++) {
      if(i==0) {
        printf("First halt for %d\n", my_rank);
        //sleep(2);
      } else if (i == 50) {
        printf("Second halt for %d\n", my_rank);
        //sleep(2);
      } else {
        printf("Printing from %d\n", my_rank);
        //sleep(1);
      }
    }

    MPI_Finalize();
    return 0;
}
