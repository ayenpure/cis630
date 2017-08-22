#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

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

void calculateRanksForRound(unordered_map<int,double> &newnodetorank,
	unordered_map<int, double> &nodetorank,
	unordered_map<int,vector<int>> &nodeadjacency) {

	unordered_map<int,vector<int>>::iterator nodeiter = nodeadjacency.begin();
	while(nodeiter != nodeadjacency.end()) {
		/*
		 * This vector has all the adjacencies.
		 */
		vector<int> neighbors = nodeiter->second;
		double rank = 0;
		for(int i = 0; i < neighbors.size(); i++) {
			int nodeId = neighbors[i];
			double neighbor_rank = nodetorank.find(nodeId)->second;
			double neighbor_degree = nodeadjacency.find(nodeId)->second.size();
			rank += neighbor_rank / neighbor_degree;
		}
		newnodetorank.insert(make_pair(nodeiter->first, rank));
		++nodeiter;
	}
}

int main(int argc, char* argv[]) {
	/*high_resolution_clock::time_point start,end;
	duration<double> required_time;
	start = high_resolution_clock::now();*/
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
	/*end = high_resolution_clock::now();
	required_time = duration_cast<duration<double>>(end - start);
	cout << "time required to read graph : " << required_time.count() << endl;
	start.~time_point();
	end.~time_point();
	required_time.~duration();*/
	initRanks(nodeadjacency, nodetorank);
	int numrounds = atoi(argv[2]);
	for(int round = 0;round < numrounds; round++) {
		unordered_map<int, double> roundnodetorank;
		//start = high_resolution_clock::now();
		calculateRanksForRound(roundnodetorank, nodetorank, nodeadjacency);
		//end = high_resolution_clock::now();
		//cout << "time required for round " << round+1 << " : "<< required_time.count() << endl;
		/*start.~time_point();
		end.~time_point();
		required_time.~duration();*/
		nodetorank = roundnodetorank;
	}
	return 0;
}
