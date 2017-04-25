#include <iostream>
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
  cout << "Returning max nodes " << max_node << endl;
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

void read_graph_data(char *file_name, int** abjacency_list, int max_node) {
  int nodeId1, nodeId2;
  int incidences[max_node+1] = {0};
  ifstream toRead(file_name);
  if (toRead.is_open()) {
  	while ( toRead >> nodeId1 >> nodeId2) {
      abjacency_list[nodeId1][incidences[nodeId1]] = nodeId2;
      abjacency_list[nodeId2][incidences[nodeId2]] = nodeId1;
      incidences[nodeId1]++;
      incidences[nodeId2]++;
  	}
  	toRead.close();
  }
}

int** allocatememory(int node_degrees[], int max_node) {
  int** adjacency_lists;
  adjacency_lists = (int**)malloc((max_node+1)*(sizeof(int*)));
  for(int i = 1; i <= max_node; i++) {
    adjacency_lists[i] = (int*)malloc(node_degrees[i]*sizeof(int));
  }
  return adjacency_lists;
}

int main(int argc,char* argv[]) {
  int number_of_rounds = atoi(argv[2]);
  const int max_node = get_max_node(argv[1]);
  int node_degrees[max_node+1]= {0};
  get_node_degrees(argv[1], node_degrees);
  /*for(int i = 1; i <= max_node; i++) {
    cout << "Degree for " << i << " is "<< node_degrees[i] << endl;
  }*/
  int** abjacency_list = allocatememory(node_degrees, max_node);
  read_graph_data(argv[1],abjacency_list, max_node);
  cout << "memory allocated successfully" << endl;
  return 0;
}
