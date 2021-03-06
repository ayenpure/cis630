#include "mpi.h"
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <ratio>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/io.h>
#include <sys/mman.h>
#define  MASTER	0

using namespace std;
using namespace std::chrono;

#define REDUCE 100
#define SCATTER 200

void getMaxNodeAndEdges(char* nodeInfoFile, int *maxNode, int* maxEdges) {
  *maxNode = 0;
  *maxEdges = 0;
  int i,snode, sdegree, spartition;
  /*ifstream toRead(nodeInfoFile);
  if (toRead.is_open()) {
  	while (toRead >> snode >> sdegree >> spartition) {
      if(snode > *maxNode)
        *maxNode = snode;
      *maxEdges+=sdegree;
  	}
  	toRead.close();
  }*/
  struct stat s;
  int fd = open (nodeInfoFile, O_RDONLY);
  int status = fstat (fd, & s);
  /* Get the size of the file. */
  int size = s.st_size;
  char *file = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
  char *lstart = file, *lend;
  for (i = 0; i < size;i++) {
    if(file[i] == '\n') {
      snode = strtol(lstart, &lend, 10);
      sdegree = strtol(lend, &lend, 10);
      spartition = strtol(lend, NULL, 10);
      lstart = file + i + 1;
      if(snode > *maxNode)
        *maxNode = snode;
      *maxEdges += sdegree;
    }
  }
  munmap(file, size);
  *maxEdges /= 2;
}

void writeToFile(double* roundRanks, int numberOfRounds, int* nodeDegree, int* isLocalNode, int allocationSize, int rank) {
    string name = to_string(rank) + ".txt";
    int i,j;
    FILE* fout = fopen(name.c_str(), "w");
    for (i = 1; i < allocationSize; i++)
    {
      if(isLocalNode[i] == rank){
        fprintf(fout, "%d\t%d\t", i, nodeDegree[i]);
        for(int j = 1; j <= numberOfRounds; j++) {
          fprintf(fout, "%lf\t", roundRanks[(j*allocationSize)+i]);
        }
        fprintf(fout, "\n");
      }
    }
    fclose(fout);
}

void getNodeInfo(char *nodeInfoFile, int* nodeDegree, int* isLocalNode, int rank) {
  int i,snode, sdegree, srank;
  /*ifstream nodeInfo(nodeInfoFile);
  if(nodeInfo.is_open()) {
    while(nodeInfo >> snode >> sdegree >> srank) {
        isLocalNode[snode] = srank;
        nodeDegree[snode] = sdegree;
    }
    nodeInfo.close();
  }*/
  struct stat s;
  int fd = open (nodeInfoFile, O_RDONLY);
  int status = fstat (fd, & s);
  /* Get the size of the file. */
  int size = s.st_size;
  char *file = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
  char *lstart = file, *lend;
  for (i = 0; i < size;i++) {
    if(file[i] == '\n') {
      snode = strtol(lstart, &lend, 10);
      sdegree = strtol(lend, &lend, 10);
      srank = strtol(lend, NULL, 10);
      lstart = file + i + 1;
      isLocalNode[snode] = srank;
      nodeDegree[snode] = sdegree;
    }
  }
  munmap(file, size);
}

void getEdgeInfo(char *edgeListFile, int** edgeList) {
  int snode, dnode, edge = 0;
  /*ifstream edgeInfo(edgeListFile);
  if(edgeInfo.is_open()) {
    while(edgeInfo >> snode >> dnode) {
      edgeList[i][0] = snode;
      edgeList[i][1] = dnode;
      i++;
    }
    edgeInfo.close();
  }*/
  struct stat s;
  int fd = open (edgeListFile, O_RDONLY);
  int status = fstat (fd, & s);
  /* Get the size of the file. */
  int size = s.st_size;
  char *file = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
  char *lstart = file, *lend;
  for (int i = 0; i < size;i++) {
    if(file[i] == '\n') {
      snode = strtol(lstart, &lend, 10);
      dnode = strtol(lend, NULL, 10);
      lstart = file + i + 1;
      edgeList[edge][0] = snode;
      edgeList[edge][1] = dnode;
      edge++;
    }
  }
  munmap(file, size);
}

void printTime(high_resolution_clock::time_point start,
  high_resolution_clock::time_point end, string message, string rmessage,
  int rank, int reduce, int reduceonly, int numProcesses) {
  double actDuration;
  duration<double> tspan;
  tspan = duration_cast<duration<double>>(end - start);
  actDuration = tspan.count();
  if(!reduceonly)
    cout << message << " " << rank << " : " << actDuration << " sec." << endl;
  if(reduce) {
    if(rank == MASTER) {
      //MPI_Recv
      double maxTime = actDuration, toReceive = 0;
      for(int i = 1; i < numProcesses; i++) {
        MPI_Recv(&toReceive, 1, MPI_DOUBLE, i , REDUCE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(toReceive > maxTime) {
          maxTime = toReceive;
        }
      }
      cout << rmessage << " : " << maxTime << " sec." << endl;
    } else {
      MPI_Send(&actDuration, 1, MPI_DOUBLE, MASTER, REDUCE, MPI_COMM_WORLD);
    }
  }
}

void EXCHANGE(double* reduce, int LIMIT, int* nodeLocation, int rank, int roundid, int numProcesses) {
  int i,j;
  int scattertag = roundid + SCATTER;
  int reducetag = roundid + REDUCE;
  double* collect = new double[LIMIT];
  //Reduce
  if(rank == MASTER) {
    //only receive
    for(i = 1; i < numProcesses; i++) {
      MPI_Recv(&collect[0], LIMIT, MPI_DOUBLE, i, reducetag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(j = 0; j < LIMIT; j++) {
        reduce[j] += collect[j];
      }
    }
  } else {
    //only send
    MPI_Send(&reduce[0], LIMIT, MPI_DOUBLE, MASTER, reducetag, MPI_COMM_WORLD);
  }
  //Scatter
  if(rank == MASTER) {
    for(i = 1; i < numProcesses; i++) {
      MPI_Send(&reduce[0], LIMIT, MPI_DOUBLE, i, scattertag, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(&reduce[0], LIMIT, MPI_DOUBLE, MASTER, scattertag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

int main (int argc, char *argv[]) {
  high_resolution_clock::time_point start, end, tstart, tend;
  tstart = high_resolution_clock::now();
  char* nodeInfoFile = argv[1];
  char* edgeListFile = argv[2];

  int  numProcesses, rank;
  int i, j, snode, dnode;
  string message, rmessage;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

  start = high_resolution_clock::now();
  const int numberOfRounds = atoi(argv[3]);
  int maxNode, maxEdges;
  getMaxNodeAndEdges(nodeInfoFile, &maxNode, &maxEdges);
  const int allocationSize = maxNode + 1;
  int *nodeDegree = new int[allocationSize];
  int *isLocalNode = new int[allocationSize];
  getNodeInfo(nodeInfoFile, nodeDegree, isLocalNode, rank);
  int **edgeList = new int*[maxEdges];
  for(i = 0; i < maxEdges; i++)
    edgeList[i] = new int[2];
  getEdgeInfo(edgeListFile, edgeList);
  end = high_resolution_clock::now();
  message = "-- Time to read input files, partition";
  printTime(start, end, message, rmessage, rank, 0, 0, numProcesses);

  double* reduce = new double[allocationSize];
  double* roundRanks;
  roundRanks = new double[(numberOfRounds + 1)*allocationSize];
  for (i = 0; i < allocationSize; i++) {
    roundRanks[i] = 1;
  }

  for(i = 1; i <= numberOfRounds; i++) {
    start = high_resolution_clock::now();
    for (j = 0; j <= maxNode; j++) {
      reduce[j] = 0;
    }
    for(j = 0; j < maxEdges; j++) {
      snode = edgeList[j][0];
      if(isLocalNode[snode] == rank) {
        dnode = edgeList[j][1];
        reduce[snode] += roundRanks[(i-1)*allocationSize + dnode] / nodeDegree[dnode];
        reduce[dnode] += roundRanks[(i-1)*allocationSize + snode] / nodeDegree[snode];
      }
    }
    end = high_resolution_clock::now();
    EXCHANGE(reduce, allocationSize, isLocalNode, rank, i, numProcesses);
    copy(reduce, reduce + allocationSize, &roundRanks[i*allocationSize]);
    message = "-- Time for round " + to_string(i) +  ", partition";
    rmessage = "Total time for round " + to_string(i);
    printTime(start, end, message, rmessage, rank, 1, 0, numProcesses);
  }
  start = high_resolution_clock::now();
  writeToFile(roundRanks, numberOfRounds,nodeDegree, isLocalNode, allocationSize, rank);
  end = high_resolution_clock::now();
  message = "-- Time to write output file, partition";
  rmessage = "Total time to write output files ";
  printTime(start, end, message, rmessage, rank, 1, 0, numProcesses);
  tend = high_resolution_clock::now();
  rmessage = "Total time for execution of program ";
  printTime(tstart, tend, message, rmessage, rank, 1, 1, numProcesses);
  MPI_Finalize();
}
