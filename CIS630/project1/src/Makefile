all: pagerank mpi fileread sr
pagerank : pagerankflatarray.cxx
	g++ -g -std=c++11 -O3 -o pagerank pagerankflatarray.cxx
mpi : pagerankmpi.cxx
	mpicxx -g -std=c++11 -O3 -o mpipagerank pagerankmpi.cxx
	#g++ -g -std=c++11 -O0 -o mpipagerank pagerankmpi.cxx
sr : pagerankmpisr.cxx
	mpicxx -g -std=c++11 -O3 -o srmpipagerank pagerankmpisr.cxx
fileread : fileread.cxx
	g++ -g -std=c++11 -O0 -o fileread fileread.cxx
