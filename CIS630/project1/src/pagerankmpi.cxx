/******************************************************************************
* FILE: mpi_helloBsend.c
* DESCRIPTION:
*   MPI tutorial example code: Simple hello world program that uses blocking
*   send/receive routines.
* AUTHOR: Blaise Barney
* LAST REVISED: 06/08/15
******************************************************************************/
#include "mpi.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#define  MASTER	0

using namespace std;

int getMaxNode(char* nodeInfoFile) {
  int maxNode = 0, nodeId1, nodeId2;
  ifstream toRead(nodeInfoFile);
  if (toRead.is_open()) {
  	while (toRead >> nodeId1 >> nodeId2) {
      if(nodeId1 > maxNode)
        maxNode = nodeId1;
      if(nodeId2 > maxNode)
        maxNode = nodeId2;
  	}
  	toRead.close();
  }
  return maxNode;
}

void writeToFile(double** roundRanks, int numberOfRounds, int* nodeDegree, int* isLocalNode, int maxNode, int rank) {
    int i;
    ofstream toWrite;
    string name = to_string(rank) + ".txt";
		toWrite.open(name);
    toWrite << fixed;
    toWrite << setprecision(6);
		for(i = 1; i <= maxNode; i++) {
      if(isLocalNode[i]){
        toWrite << i << "\t" << nodeDegree[i] << "\t";
        for(int j = 1; j <= numberOfRounds; j++) {
            toWrite << roundRanks[j][i] << "\t";
        }
        toWrite << endl;
      }
    }
		toWrite.close();
}

void getNodeInfo(char *nodeInfoFile, int* nodeDegree, int* nodeLocation, int** adjacencyList, int rank) {
  int snode, sdegree, srank;
  ifstream nodeInfo(nodeInfoFile);
  if(nodeInfo.is_open()) {
    while(nodeInfo >> snode >> sdegree >> srank) {
      if(srank != rank) {
          nodeLocation[snode] = 0;
          nodeDegree[snode] = sdegree;
          continue;
      }
      nodeLocation[snode] = 1;
      nodeDegree[snode] = sdegree;
      adjacencyList[snode] = new int[sdegree];
    }
    nodeInfo.close();
  }
}

void addEdgeBetweenNodes(int snode, int dnode,  int** adjacencyList, int index) {
  adjacencyList[snode][index] = dnode;
}

void getEdgeInfo(char *edgeListFile, int** adjacencyList, int* isLocalNode, int allocationSize) {
  int snode, dnode;
  int i;
  int index[allocationSize];
  for(i = 0; i < allocationSize; i++) {
    index[i] = 0;
  }
  ifstream edgeList(edgeListFile);
  if(edgeList.is_open()) {
    while(edgeList >> snode >> dnode) {
      if(isLocalNode[snode]){
        addEdgeBetweenNodes(snode, dnode, adjacencyList, index[snode]);
        ++index[snode];
      }
      if(isLocalNode[dnode]) {
        addEdgeBetweenNodes(dnode, snode, adjacencyList, index[dnode]);
        ++index[dnode];
      }
    }
    edgeList.close();
  }
}

int main (int argc, char *argv[])
{
  int  numProcesses, rank;
  int i, j, snode, dnode;
  double sum;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

  //read file and populate data structures for the rank
  char* nodeInfoFile = argv[1];
  char* edgeListFile = argv[2];
  const int numberOfRounds = atoi(argv[3]);
  const int maxNode = getMaxNode(nodeInfoFile);
  const int allocationSize = maxNode + 1;
  int *nodeDegree = new int[allocationSize];
  int *isLocalNode = new int[allocationSize];;
  int **adjacencyList = new int*[allocationSize];
  getNodeInfo(nodeInfoFile, nodeDegree, isLocalNode, adjacencyList, rank);
  getEdgeInfo(edgeListFile, adjacencyList, isLocalNode, allocationSize);

  //allocate space for rounds of page credits calculation
  //either a 2D array or a map of vectors storing the credits
  double** roundRanks;
  double* reduce = new double[allocationSize];
  roundRanks = new double*[numberOfRounds + 1];
  for (i = 0; i <= numberOfRounds; i++) {
    roundRanks[i] = new double[maxNode];
  }
  for (i = 0; i <= maxNode; i++) {
    roundRanks[0][i] = 1;
  }
  /*start rounds {
    loop through nodes of this node.
    If the other node is on this machine get credits and degree.
    If the other node is not on this machine, request its credits
    MPI_Send()/MPI_Rec()
    -the tag be the round ID
    -message req be an int and receive as an int
    calculate the rank, store the result, move on.
  }*/
  for(i = 1; i <= numberOfRounds; i++) {
    for (j = 0; j <= maxNode; j++) {
      reduce[j] = 0;
    }
    //go to every local node
    for(snode = 1; snode <= maxNode; snode++) {
      if(!isLocalNode[snode])
        continue;
      sum = 0;
      //for each neighbor locally available calculate locally the credits.
      for(j = 0; j < nodeDegree[snode]; j++) {
        dnode = adjacencyList[snode][j];
        /*if(!isLocalNode[dnode])
          continue;*/
        sum += roundRanks[i-1][dnode] / nodeDegree[dnode];
      }
      reduce[snode] = sum;
    }
    //store these credits in the designates arrays
    MPI_Allreduce(reduce, roundRanks[i], allocationSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /*if(rank == MASTER)*/
  }
  writeToFile(roundRanks, numberOfRounds,nodeDegree, isLocalNode, maxNode, rank);
  /*
    Finish up the rounds, collect all the credits in the master
    write them to the output file
  */
  MPI_Finalize();
}
