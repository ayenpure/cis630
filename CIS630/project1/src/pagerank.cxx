#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <chrono>
#include <map>
#include <unordered_map>
//1861227
//node232957
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
	cout << "initializing ranks" << endl;
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
		cout << "Calculating new page ranks" << endl;
		unordered_map<int,vector<int>>::iterator nodeiter = nodeadjacency.begin();
		while(nodeiter != nodeadjacency.end()) {
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

void writeToFile(unordered_map<int, double> &nodetorank,
	unordered_map<int,vector<int>> &nodeadjacency) {
		cout << "writing to file" << endl;
		ofstream toWrite;
		toWrite.open("output.txt");
		std::map<int,vector<int>> ordered(nodeadjacency.begin(), nodeadjacency.end());
		auto nodeiter = ordered.begin();
		while(nodeiter != ordered.end()) {
			toWrite << nodeiter->first << "\t" << nodeiter->second.size() << "\t" << nodetorank.find(nodeiter->first)->second << endl;
			++nodeiter;
		}
		cout << nodetorank.size() << endl;
		toWrite.close();
}

int main(int argc, char* argv[]) {
	/*high_resolution_clock::time_point start,end;
	duration<double> required_time;
	start = high_resolution_clock::now();*/
	string line;
  ifstream toRead(argv[1]);
	unordered_map<int, double> nodetorank;
	unordered_map<int,vector<int>> nodeadjacency;
	int nodeId1,nodeId2;
 	if (toRead.is_open()) {
		while ( toRead >> nodeId1 >> nodeId2) {
			addEdgeBetweenNodes(nodeId1, nodeId2, nodeadjacency);
			addEdgeBetweenNodes(nodeId2, nodeId1, nodeadjacency);
		}
		toRead.close();
	}
	/*end = high_resolution_clock::now();
	required_time = duration_cast<duration<double>>(end - start);
	cout << "time required to read graph : " << required_time.count() << endl;
	start.~time_point();
	end.~time_point();
	required_time.~duration();*/
	initRanks(nodeadjacency, nodetorank);
	int numrounds = atoi(argv[2]);
	unordered_map<int, double> roundnodetorank;
	for(int round = 0;round < numrounds; round++) {
		//start = high_resolution_clock::now();
		calculateRanksForRound(roundnodetorank, nodetorank, nodeadjacency);
		//end = high_resolution_clock::now();
		//cout << "time required for round " << round+1 << " : "<< required_time.count() << endl;
		/*start.~time_point();
		end.~time_point();
		required_time.~duration();*/
	}
	writeToFile(roundnodetorank,nodeadjacency);
}
