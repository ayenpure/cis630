#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>

using namespace::std;

void addEdgeBetweenNodes(int nodeId1, int nodeId2, unordered_map<int,vector<int>> &nodeadjacency) {
	std::unordered_map<int,vector<int>>::iterator nodeiter;
	nodeiter = nodeadjacency.find(nodeId1);
	if(nodeiter != nodeadjacency.end()) {
		nodeiter->second.push_back(nodeId2);
	} else {
		vector<int> adjlist;
		adjlist.push_back(nodeId2);
		nodeadjacency.insert(make_pair(nodeId1, adjlist));
	}
}

int main(int argc, char* argv[]) {
	string line;
  ifstream graph(argv[1]);
	unordered_map<int, int> nodetorank;
	unordered_map<int,vector<int>> nodeadjacency;
	int nodeId1,nodeId2;
	//stringstream toRead;
 	if (graph.is_open()) {
		while ( /*getline (graph,line)*/
	 					graph >> nodeId1 >> nodeId2) {
						/*toRead.str(line);
						toRead.clear();
						toRead >> nodeId1 >> nodeId2;*/
						//cout << nodeId1 << " " << nodeId2 << endl;
						addEdgeBetweenNodes(nodeId1, nodeId2, nodeadjacency);
						addEdgeBetweenNodes(nodeId2, nodeId1, nodeadjacency);
		}
		graph.close();
	}
	/*std::unordered_map<int,vector<int>>::iterator nodeiter;
	for() {

	}*/
	return 0;
}
