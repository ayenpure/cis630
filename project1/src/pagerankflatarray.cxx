#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

int get_max_node(char* file_name) {
  int max_node = 0, nodeId1, nodeId2;
  ifstream toRead(file_name);
  if (toRead.is_open()) {
  	while (toRead >> nodeId1 >> nodeId2) {
      if(nodeId1 > max_node)
        max_node = nodeId1;
      else if(nodeId2 > max_node)
        max_node = nodeId2;
  	}
  	toRead.close();
  }
  //cout << "Returning max nodes " << max_node << endl;
  return max_node;
}

void get_node_degrees(char *file_name, int node_degrees[]) {
  int nodeId1, nodeId2;
  ifstream toRead(file_name);
  if (toRead.is_open()) {
  	while ( toRead >> nodeId1 >> nodeId2) {
      node_degrees[nodeId1]++;
      node_degrees[nodeId2]++;
  	}
  	toRead.close();
  }
}

void read_graph_data(char *file_name, int** adjacency_list, int allocation_size) {
  int nodeId1, nodeId2;
  int *incidences = new int[allocation_size];
  ifstream toRead(file_name);
  if (toRead.is_open()) {
  	while ( toRead >> nodeId1 >> nodeId2) {
      adjacency_list[nodeId1][incidences[nodeId1]] = nodeId2;
      adjacency_list[nodeId2][incidences[nodeId2]] = nodeId1;
      incidences[nodeId1]++;
      incidences[nodeId2]++;
  	}
  	toRead.close();
  }
  delete [] incidences;
}

int** allocate_memory_for_adjacency(int node_degrees[], int allocation_size) {
  int** adjacency_lists;
  adjacency_lists = (int**)malloc(allocation_size*sizeof(int*));
  for(int i = 1; i <= allocation_size; i++) {
    adjacency_lists[i] = (int*)malloc(node_degrees[i]*sizeof(int));
  }
  return adjacency_lists;
}

double** allocate_memory_for_rounds(int number_of_rounds, int allocation_size) {
  double** round_ranks;
  round_ranks = (double**)malloc((number_of_rounds+1)*sizeof(double*));
  for(int i = 0; i <= number_of_rounds; i++) {
    round_ranks[i] = (double*)malloc(allocation_size*sizeof(double));
  }
  return round_ranks;
}

void initRanks(double* round_rank, int allocation_size) {
  for(int i = 1; i<=allocation_size; i++) {
    round_rank[i] = 1;
  }
}

void performPageRankRounds(double** round_ranks, int** adjacency_list, int* node_degrees, int allocation_size, int number_of_rounds) {
  for(int i = 1; i <= number_of_rounds; i++) {
    for(int j = 1; j <= allocation_size; j++) {
      double curr_rank = 0;
      for(int k = 0; k < node_degrees[j]; k++) {
        int curr_inc = adjacency_list[j][k];
        curr_rank += round_ranks[i-1][curr_inc] / node_degrees[curr_inc];
      }
      round_ranks[i][j] = curr_rank;
    }
  }
}

void writeToFile(int *node_degrees, double** round_ranks,int allocation_size, int number_of_rounds) {
  ofstream toWrite;
  toWrite.open("output.txt");
  toWrite << std::setprecision(6);
  toWrite << std::fixed;
  for(int i = 1; i <= allocation_size; i++) {
    if(node_degrees[i] == 0)
      continue;
    toWrite << i << "\t" << node_degrees[i] << "\t";
    for (int j = 1; j <= number_of_rounds; j++) {
      toWrite << round_ranks[j][i] << "\t";
    }
    toWrite << endl;
  }
  toWrite.close();
}

int main(int argc,char* argv[]) {
  if(argc != 3) {
    cerr << "Usage : <pagerank> <input graph> <number of rounds>" << endl;
    exit(EXIT_FAILURE);
  }
  int number_of_rounds = atoi(argv[2]);
  const int max_node = get_max_node(argv[1]);
  const int allocation_size = max_node + 1;
  int *node_degrees = new int[allocation_size];
  get_node_degrees(argv[1], node_degrees);
  int** adjacency_list = allocate_memory_for_adjacency(node_degrees, allocation_size);
  read_graph_data(argv[1],adjacency_list, allocation_size);
  double** round_ranks =  allocate_memory_for_rounds(number_of_rounds+1, allocation_size);
  initRanks(round_ranks[0],allocation_size);
  performPageRankRounds(round_ranks, adjacency_list, node_degrees, allocation_size, number_of_rounds);
  writeToFile(node_degrees, round_ranks, allocation_size, number_of_rounds);
  return 0;
}
