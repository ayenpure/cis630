all: mpi comment
mpi : pagerankmpi.cxx
	@mpicxx -g -std=c++11 -O3 -o mpipagerank pagerankmpi.cxx
comment :
	@echo ________________________________________________________________________________________
	@echo
	@echo Usage : mpipagerank \<graph file\> \<partition file\> \<num rounds\> \<num partitions\>
	@echo ________________________________________________________________________________________
run :
	mpiexec -np 2 ./mpipagerank /cs/classes/www/17S/cis630/LectureNotes/fl_compact.tab /cs/classes/www/17S/cis630/LectureNotes/fl_compact_part.2 5 2
clean : 
	@rm -f mpipagerank
