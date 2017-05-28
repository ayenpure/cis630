/*
 *  Abhishek Yenpure
 *  CIS 630
 *  Project 1
 */
#include <chrono>
#include <ctime>
#include <iostream>
#include <mpi.h>
#include <ratio>
#include <stdio.h>
#include <stdlib.h>
#define  MASTER	0

using namespace std;
using namespace std::chrono;

/* 
I have used this function as a helper and have used it from a stackoverflow solution
http://stackoverflow.com/questions/16826422/c-most-efficient-way-to-convert-string-to-int-faster-than-atoi
*/
int fast_atoi(const char * str) {
    int val = 0;
    while( *str ) {
        val = val*10 + (*str++ - '0');
    }
    return val;
}

/*
This function does the job of fetching the highest node and the number of edges of the graph
*/
void getMaxNodeAndEdges(const char* nodeInfoFile, int *maxNode, int* maxEdges) {
  *maxNode = 0;
  *maxEdges = 0;
  size_t readBytes;
  char readNumber[11];
  int i, nodeInfo[3], readlength = 0, numCount = 0;
  char *buffer = (char*)malloc(1024*1024);
  FILE *inputFile = fopen(nodeInfoFile, "r");
  if(inputFile == NULL) {
    cerr << "Error reading the input file : " << nodeInfoFile << endl;
    exit(1);
  }
  while(!feof(inputFile)) {
    readBytes = fread(buffer, 1, 1024*1024, inputFile);
    for(i = 0; i < readBytes; i++) {
      if(isdigit(buffer[i])) {
        readNumber[readlength++] = buffer[i];
      } else if(buffer[i] == '\n' || buffer[i] == '\t' || buffer[i] == EOF ) {
        readNumber[readlength] = '\0';
        readlength = 0;
        nodeInfo[numCount++] = fast_atoi(readNumber);
        if(numCount == 3) {
          numCount = 0;
          if(nodeInfo[0] > *maxNode)
            *maxNode = nodeInfo[0];
          *maxEdges += nodeInfo[1];
        }
      }
    }
  }
  free(buffer);
  fclose(inputFile);
  *maxEdges /= 2;
}

/*
This function updates the data strucutres storing node locations and degrees
*/
void getNodeInfo(const char *nodeInfoFile, int* nodeDegree, int* isLocalNode, int rank) {
  size_t readBytes;
  char readNumber[11];
  int i, nodeInfo[3], readlength = 0, numCount = 0;
  char *buffer = (char*)malloc(1024*1024);
  FILE *inputFile = fopen(nodeInfoFile, "r");
  if(inputFile == NULL) {
    cerr << "Error reading the input file : " << nodeInfoFile << endl;
    exit(1);
  }
  while(!feof(inputFile)) {
    readBytes = fread(buffer, 1, 1024*1024, inputFile);
    for(i = 0; i < readBytes; i++) {
      if(isdigit(buffer[i])) {
        readNumber[readlength++] = buffer[i];
      } else if(buffer[i] == '\n' || buffer[i] == '\t' || buffer[i] == EOF ) {
        readNumber[readlength] = '\0';
        readlength = 0;
        nodeInfo[numCount++] = fast_atoi(readNumber);
        if(numCount == 3) {
          numCount = 0;
          isLocalNode[nodeInfo[0]] = nodeInfo[2];
          nodeDegree[nodeInfo[0]] = nodeInfo[1];
        }
      }
    }
  }
  free(buffer);
  fclose(inputFile);
}

/*
This function fetches the edgelist from the graph
*/
void getEdgeInfo(const char *edgeListFile, int** edgeList) {
  size_t readBytes;
  char readNumber[11];
  int i, nodeInfo[2], readlength = 0, numCount = 0, edge = 0;
  char *buffer = (char*)malloc(1024*1024);
  FILE *inputFile = fopen(edgeListFile, "r");
  if(inputFile == NULL) {
    cerr << "Error reading the input file : " << edgeListFile << endl;
    exit(1);
  }
  while(!feof(inputFile)) {
    readBytes = fread(buffer, 1, 1024*1024, inputFile);
    for(i = 0; i < readBytes; i++) {
      if(isdigit(buffer[i])) {
        readNumber[readlength++] = buffer[i];
      } else if(buffer[i] == '\n' || buffer[i] == '\t' || buffer[i] == EOF ) {
        readNumber[readlength] = '\0';
        readlength = 0;
        nodeInfo[numCount++] = fast_atoi(readNumber);
        if(numCount == 2) {
          numCount = 0;
          edgeList[edge][0] = nodeInfo[0];
          edgeList[edge][1] = nodeInfo[1];
          edge++;
        }
      }
    }
  }
  free(buffer);
  fclose(inputFile);
}

/*
This function writes the output of all rounds to an output file
*/
void writeToFile(double* roundRanks, int numberOfRounds, int* nodeDegree, int* isLocalNode, int allocationSize, int rank) {
  string name = to_string(rank) + ".out";
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

/*
This function prints the times to the output
*/
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
    double maxTime;
    MPI_Reduce(&actDuration, &maxTime, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    if(rank == MASTER)
      cout << rmessage << " : " << maxTime << " sec." << endl;
  }
}


int main (int argc, char *argv[]) {
  if(argc != 5) {
    cerr << "Err: Usage : " << argv[0] << " <graph file> <partition file> <number of rounds> <number of partitions>" << endl;
    exit(1);
  }
  high_resolution_clock::time_point start, end, tstart, tend;
  tstart = high_resolution_clock::now();
  const char* edgeListFile = argv[1];
  const char* nodeInfoFile = argv[2];
  const int numberOfRounds = fast_atoi(argv[3]);
  const int numPartitions = fast_atoi(argv[4]);
  int  numProcesses, rank;
  int i, j, snode, dnode;
  string message, rmessage;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  if ( numProcesses != numPartitions) {
    cerr << "Err : Number of partitions was not equal to the number of spawned MPI processes" << endl;
    exit(1);
  }
  start = high_resolution_clock::now();
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

  /*
	Begin rounds
  */
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
    MPI_Allreduce(reduce, &roundRanks[i*allocationSize], allocationSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
