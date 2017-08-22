all: pagerank mpi mpifread
pagerank : pagerankflatarray.cxx
	g++ -g -std=c++11 -O3 -o pagerank pagerankflatarray.cxx
mpi : pagerankmpi.cxx
	mpicxx -g -std=c++11 -O3 -o mpipagerank pagerankmpi.cxx
	#g++ -g -std=c++11 -O0 -o mpipagerank pagerankmpi.cxx

mpifread : pagerankmpiparse.cxx
	mpicxx -g -std=c++11 -O3 -o freadpagerank pagerankmpiparse.cxx
	#g++ -g -std=c++11 -O0 -o mpipagerank pagerankmpi.cxx
