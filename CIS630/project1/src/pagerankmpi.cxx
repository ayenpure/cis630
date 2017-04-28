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
    /*snodeDegrees = new int[maxNode+1];
    for(int i = 0; i <= maxNode; i++) {
      snodeDegrees[i] = 0;
    }*/
    //int indextoinsert = 0;
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
  unordered_map<int, int> nodeLocation;
  unordered_map<int, int> nodeDegree;
  unordered_map<int, vector<int>> nodeAdjacencies;
  int maxNode = readAndPopulate(nodeInfoFile, edgeListFile, rank, nodeLocation, nodeDegree, nodeAdjacencies);
  cout << nodeLocation.size() << " " << nodeDegree.size() << " " << nodeAdjacencies.size() << endl;

  if(rank == MASTER) {
    map<int,vector<int>> ordered(nodeAdjacencies.begin(), nodeAdjacencies.end());
    ofstream toWrite("verify");
    toWrite << fixed;
    for(auto iter = ordered.begin(); iter != ordered.end(); iter++) {
      toWrite << iter->first << "\t" << (nodeDegree.find(iter->first))->second << " : ";
      for(int i = 0; i < iter->second.size(); i++) {
        toWrite << iter->second[i] << "\t";
      }
      toWrite << endl;
    }
    toWrite.close();
  }

  //allocate space for rounds of page credits calculation
  //either a 2D array or a map of vectors storing the credits
  int** roundRanks = new int*[numberOfRounds + 1];
  int* remoteRanks = new int[maxNode +1];
  for(int i = 0; i <= numberOfRounds; i++) {
    roundRanks = new int[maxNode];
  }
  for (int i = 0; i <= maxNode; i++) {
    roundRanks[0][i] = 1;
    remoteRanks[i] = -1;
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

  }
  /*
    Finish up the rounds, collect all the credits in the master
    write them to the output file
  */

  MPI_Finalize();
}
