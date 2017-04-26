all: pagerank mpi
pagerank : pagerankflatarray.cxx
	g++ -g -std=c++11 -O3 -o pagerank pagerankflatarray.cxx

mpi : pagerankmpi.cxx
	mpicxx -g -std=c++11 -O3 -o mpipagerank pagerankmpi.cxx
