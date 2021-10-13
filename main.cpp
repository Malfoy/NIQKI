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
    // monidex.insert_file_lines("test.fa");
    // monidex.query_file_lines("test.fa");
    // monidex.insert_file_lines("1000Bact.fa");
    monidex.insert_file_of_file_whole("fof1000.txt");
    monidex.query_file_of_file_whole("fof1000.txt");
    
    // monidex.query_file_lines("random_1mb_1k.fa");
    // monidex.insert_file_lines("random_1mb_10.fa");
    // monidex.query_file_lines("random_1mb_10.fa");
    return 0;
}