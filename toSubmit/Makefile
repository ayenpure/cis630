all: mpi comment
mpi : pagerankmpi.cxx
	mpicxx -g -std=c++11 -O3 -o mpipagerank pagerankmpi.cxx
comment :
	@echo ________________________________________________________________________________________
	@echo
	@echo Usage : mpipagerank \<graph file\> \<partition file\> \<num rounds\> \<num partitions\>
	@echo ________________________________________________________________________________________
run :
	mpiexec -np 2 ./mpipagerank fl_compact.tab fl_compact_part.2 5 2
