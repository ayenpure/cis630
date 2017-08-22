#include <iostream>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/io.h>
#include <sys/mman.h>
#include <chrono>
#include <ctime>
#include <ratio>
#include <stdlib.h>
#include <fstream>

using namespace std;
using namespace std::chrono;

void printTime(high_resolution_clock::time_point start,
  high_resolution_clock::time_point end) {
  duration<double> tspan;
  tspan = duration_cast<duration<double>>(end - start);
  cout << "Time required for file reading : " << tspan.count() << " sec." << endl;
}


int main(int argc, char const *argv[])
{
    high_resolution_clock::time_point start1, end1, start2, end2;
    start1 = high_resolution_clock::now();
    char *f;
    int size;
    struct stat s;
    const char * file_name = argv[1];
    int fd = open (argv[1], O_RDONLY);

    /* Get the size of the file. */
    int status = fstat (fd, & s);
    size = s.st_size;
    char intlen[10];
    int adv = 0;
    f = (char *) mmap (0, size, PROT_READ, MAP_PRIVATE, fd, 0);
    int nums[3];
    int cnt,line = 0;
    int found = 0;
    for (int i = 0; i < size;i++) {
        while(f[i] != '\t' && f[i] != '\n' && f[i] != '\0' && i < size) {
          found = 1;
          intlen[adv++] = f[i];
          i++;
        }
        intlen[adv] = '\0';
        adv = 0;
        if(found) {
          nums[cnt++] = atoi(intlen);
          if(cnt == 2) {
            cnt = 0;
            line++;
            // cout << nums[0] << "\t" << nums[1] << "\t" << nums[2];
          }
        }
    }
    end1 = high_resolution_clock::now();
    start2 = high_resolution_clock::now();
    int snode, sdegree, spartition;
    ifstream toRead(argv[1]);
    if (toRead.is_open()) {
    	while (toRead >> snode >> sdegree) {
        //cout << snode << "\t" << sdegree << "\t" << spartition;
    	}
    	toRead.close();
    }
    end2 = high_resolution_clock::now();
    printTime(start1, end1);
    printTime(start2, end2);
    return 0;
}
