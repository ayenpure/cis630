#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>

using namespace::std;

void addEdgeBetweenNodes(int nodeId1, int nodeId2,
	unordered_map<int,vector<int>> &nodeadjacency) {
	unordered_map<int,vector<int>>::iterator nodeiter;
	nodeiter = nodeadjacency.find(nodeId1);
	if(nodeiter != nodeadjacency.end()) {
		nodeiter->second.push_back(nodeId2);
	} else {
		vector<int> adjlist;
		adjlist.push_back(nodeId2);
		nodeadjacency.insert(make_pair(nodeId1, adjlist));
	}
}

void initRanks(unordered_map<int,vector<int>> &nodeadjacency,
	unordered_map<int, double> &nodetorank) {
	unordered_map<int,vector<int>>::iterator nodeiter = nodeadjacency.begin();
	while(nodeiter != nodeadjacency.end()) {
		nodetorank.insert(make_pair(nodeiter->first, 1.));
		++nodeiter;
	}
	cout << nodetorank.size() << endl;
}

void calculateRanksForRound(unordered_map<int,
	double> &nodetorank, unordered_map<int, double> &nodetorank,
	unordered_map<int,vector<int>> &nodeadjacency) {
	while(nodeiter != nodeadjacency.end()) {
		vector<int> neighbors = nodeiter->second;
		++nodeiter;
	}
}

int main(int argc, char* argv[]) {
	string line;
  ifstream graph(argv[1]);
	unordered_map<int, double> nodetorank;
	unordered_map<int,vector<int>> nodeadjacency;
	int nodeId1,nodeId2;
 	if (graph.is_open()) {
		while ( graph >> nodeId1 >> nodeId2) {
						addEdgeBetweenNodes(nodeId1, nodeId2, nodeadjacency);
						addEdgeBetweenNodes(nodeId2, nodeId1, nodeadjacency);
		}
		graph.close();
	}
	initRanks(nodeadjacency, nodetorank);
	int numrounds = atoi(argv[2]);
	for(int round = 0;round < numrounds; round++) {
		unordered_map<int, double> roundnodetorank;
		calculateRanksForRound(roundnodetorank, nodetorank, nodeadjacency);
		nodetorank = roundnodetorank;
	}
	return 0;
}
