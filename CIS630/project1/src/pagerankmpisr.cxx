#include "mpi.h"
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <ratio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#define  MASTER	0

using namespace std;
using namespace std::chrono;

#define REDUCE 100
#define SCATTER 200

void getMaxNodeAndEdges(char* nodeInfoFile, int *maxNode, int* maxEdges) {
  *maxNode = 0;
  *maxEdges = 0;
  int snode, sdegree, spartition;
  ifstream toRead(nodeInfoFile);
  if (toRead.is_open()) {
  	while (toRead >> snode >> sdegree >> spartition) {
      if(snode > *maxNode)
        *maxNode = snode;
      *maxEdges+=sdegree;
  	}
  	toRead.close();
  }
  *maxEdges /= 2;
}

void writeToFile(double** roundRanks, int numberOfRounds, int* nodeDegree, int* isLocalNode, int maxNode, int rank) {
    int i;
    ofstream toWrite;
    string name = to_string(rank) + ".txt";
		toWrite.open(name);
    toWrite << fixed;
    toWrite << setprecision(6);
		for(i = 1; i <= maxNode; i++) {
      if(isLocalNode[i] == rank){
        toWrite << i << "\t" << nodeDegree[i] << "\t";
        for(int j = 1; j <= numberOfRounds; j++) {
            toWrite << roundRanks[j][i] << "\t";
        }
        toWrite << endl;
      }
    }
		toWrite.close();
}

void getNodeInfo(char *nodeInfoFile, int* nodeDegree, int* isLocalNode, int rank) {
  int snode, sdegree, srank;
  ifstream nodeInfo(nodeInfoFile);
  if(nodeInfo.is_open()) {
    while(nodeInfo >> snode >> sdegree >> srank) {
        isLocalNode[snode] = srank;
        nodeDegree[snode] = sdegree;
    }
    nodeInfo.close();
  }
}

void getEdgeInfo(char *edgeListFile, int** edgeList) {
  int snode, dnode, i = 0;
  ifstream edgeInfo(edgeListFile);
  if(edgeInfo.is_open()) {
    while(edgeInfo >> snode >> dnode) {
      edgeList[i][0] = snode;
      edgeList[i][1] = dnode;
      i++;
    }
    edgeInfo.close();
  }
}

void printTime(high_resolution_clock::time_point start,
  high_resolution_clock::time_point end, string message, string rmessage,
  int rank, int reduce) {
  double toReduce,toCollect = 0;
  duration<double> tspan;
  tspan = duration_cast<duration<double>>(end - start);
  toReduce = tspan.count();
  cout << message << " " << rank << " : " << toReduce << " sec." << endl;
  if(reduce) {
    if(rank == MASTER) {
      //MPI_Recv
      cout << rmessage << " : " << toCollect << " sec." << endl;
    } else {
      //MPI_Send
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
  printTime(start, end, message, rmessage, rank, 0);

  double** roundRanks;
  double* reduce = new double[allocationSize];
  roundRanks = new double*[numberOfRounds + 1];
  for (i = 0; i <= numberOfRounds; i++) {
    roundRanks[i] = new double[allocationSize];
  }
  for (i = 0; i <= maxNode; i++) {
    roundRanks[0][i] = 1;
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
        reduce[snode] += roundRanks[i-1][dnode] / nodeDegree[dnode];
        reduce[dnode] += roundRanks[i-1][snode] / nodeDegree[snode];
      }
    }
    end = high_resolution_clock::now();
    EXCHANGE(reduce, allocationSize, isLocalNode, rank, i, numProcesses);
    copy(reduce, reduce + allocationSize, roundRanks[i]);
    message = "-- Time for round " + to_string(i) +  ", partition";
    rmessage = "Total time for round " + to_string(i);
    printTime(start, end, message, rmessage, rank, 1);
  }
  start = high_resolution_clock::now();
  writeToFile(roundRanks, numberOfRounds,nodeDegree, isLocalNode, maxNode, rank);
  end = high_resolution_clock::now();
  message = "-- Time to write output file, partition";
  rmessage = "Total time to write output files " + to_string(i);
  printTime(start, end, message, rmessage, rank, 1);
  MPI_Finalize();
  tend = high_resolution_clock::now();

}