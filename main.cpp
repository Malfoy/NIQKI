#include "strict_fstream.hpp"
#include "zstr.hpp"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <bitset>
#include "index.h"




using namespace std;
using namespace chrono;



int main(int argc, char ** argv){
    Index monidex(16,31,10,4);
    time_point<system_clock> start, endindex,end;
  
    start = std::chrono::system_clock::now();
    monidex.insert_file_of_file_whole("fof100.txt");
    endindex = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endindex - start;
    std::cout << "Indexing lasted " << elapsed_seconds.count() << "s\n";
    monidex.query_file_of_file_whole("fof100.txt");
    end = std::chrono::system_clock::now();
    elapsed_seconds =end-  endindex;
    cout << "Query lasted " << elapsed_seconds.count() << "s\n";
    elapsed_seconds =end-  start;
    cout<<"whole run tool took " << elapsed_seconds.count() << endl;

    return 0;
}