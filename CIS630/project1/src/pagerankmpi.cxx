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

void addEdgeBetweenNodes(int snode, int dnode, unordered_map<int, vector<int>> &nodeAdjacencies) {
  auto adjancencies = nodeAdjacencies.find(snode);
  if(adjancencies != nodeAdjacencies.end()) {
    //indextoinsert = snodeDegrees[snode];
    adjancencies->second.push_back(dnode);//[indextoinsert] = dnode;
    //snodeDegrees[snode]++;
  }
}

int readAndPopulate(char *nodeInfoFile, char* edgeListFile, int rank,
  unordered_map<int, int> &nodeLocation,
  unordered_map<int, int> &nodeDegree,
  unordered_map<int, vector<int>> &nodeAdjacencies) {
    int snode, dnode, sdegree, srank;
    //int *snodeDegrees;
    int maxNode = 0;
    /*
     Block to allocate memory for adjancencies and init the node
     locations and degrees
    */
    ifstream nodeInfo(nodeInfoFile);
    if(nodeInfo.is_open()) {
      while(nodeInfo >> snode >> sdegree >> srank) {
        nodeLocation.insert(make_pair(snode, srank));
        if(snode > maxNode)
          maxNode = snode;
        if(srank != rank)
          continue;
        nodeDegree.insert(make_pair(snode, sdegree));
        vector<int> adjancencies;
        adjancencies.reserve(sdegree);
        nodeAdjacencies.insert(make_pair(snode, adjancencies));
      }
      nodeInfo.close();
    }
    ifstream edgeList(edgeListFile);
    if(edgeList.is_open()) {
      while(edgeList >> snode >> dnode) {
        addEdgeBetweenNodes(snode, dnode, nodeAdjacencies);
        addEdgeBetweenNodes(dnode, snode, nodeAdjacencies);
      }
      edgeList.close();
    }
    return maxNode;
}

void writeToFile(double* roundRanks, int maxNode) {
		//cout << "writing to file" << endl;
		ofstream toWrite;
		toWrite.open("output.txt");
		for(int i = 1; i <= maxNode; i++) {
      toWrite << i << " : " << roundRanks[i] << endl;
    }
		//cout << nodetorank.size() << endl;
		toWrite.close();
}

int main (int argc, char *argv[])
{
  int  numProcesses, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  //read file and populate data structures for the rank
  char* nodeInfoFile = argv[1];
  char* edgeListFile = argv[2];
  int numberOfRounds = atoi(argv[3]);
  /*
    might not need this because of MPI_AllReduce
  */
  unordered_map<int, int> nodeLocation;
  unordered_map<int, int> nodeDegree;
  unordered_map<int, vector<int>> nodeAdjacencies;
  int maxNode = readAndPopulate(nodeInfoFile, edgeListFile, rank, nodeLocation, nodeDegree, nodeAdjacencies);
  cout << nodeLocation.size() << " " << nodeDegree.size() << " " << nodeAdjacencies.size() << endl;

  //allocate space for rounds of page credits calculation
  //either a 2D array or a map of vectors storing the credits
  double** roundRanks;
  double* remoteRanks;
  double* reduce = new double[maxNode +1];
  if(rank == MASTER) {
    roundRanks = new double*[numberOfRounds + 1];
    for(int i = 0; i <= numberOfRounds; i++) {
      roundRanks[i] = new double[maxNode];
    }
    for (int i = 0; i <= maxNode; i++) {
      roundRanks[0][i] = 1;
    }
  } else {
    remoteRanks = new double[maxNode +1];
    for (int i = 0; i <= maxNode; i++) {
      remoteRanks[i] = 1;
    }
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
  for(int i = 1; i <= numberOfRounds; i++) {
    for (int i = 0; i <= maxNode; i++) {
      reduce[i] = 0;
    }
    //go to every local node
    for(auto iter = nodeAdjacencies.begin() ; iter != nodeAdjacencies.end(); iter++ ) {
      int snode = iter->first;
      vector<int> incidences = iter->second;
      double sum = 0;
      //for each neighbor locally available calculate locally the credits.
      for(int j = 0; j < incidences.size(); j++) {
        int dnode = incidences[j];
        if((nodeLocation.find(dnode))->second != rank)
          continue;
        sum += ((rank == MASTER)?roundRanks[i-1][dnode]:remoteRanks[dnode])/(nodeDegree.find(dnode))->second;
      }
      //store these credits in the designates arrays
      if(rank == MASTER) {
        reduce[snode] = sum;
      } else {
        reduce[snode] = sum;
      }
    }
    //MPI_AllReduce with the data
    if(rank == MASTER) {
      MPI_Allreduce(reduce, roundRanks[i], maxNode+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } else {
      MPI_Allreduce(reduce, remoteRanks, maxNode+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    if(rank == MASTER)
      writeToFile(roundRanks[i], maxNode);
  }
  /*
    Finish up the rounds, collect all the credits in the master
    write them to the output file
  */
  MPI_Finalize();
}
