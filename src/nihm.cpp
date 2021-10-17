#include "nihm_index.h"
#include "optionparser.h"
#include "common.h"

#include <vector>
#include <string>
#include <fstream> //readline
#include <iostream> //readline
#include <csignal>
#include <algorithm>

#include <unistd.h>     // getpid(), getcwd()
#include <sys/types.h>  // type definitions, e.g., pid_t
#include <sys/stat.h>
#include <sys/wait.h>   // wait()
#include <sys/time.h>   // gettimeofday()
#include <fcntl.h>
#include <libgen.h>
#include <signal.h>     // signal name constants and kill()

#include <stdio.h>    
#include <stdint.h>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <math.h>
#include <unistd.h>
#include <bitset>

#include <zlib.h>

using namespace std;
using namespace chrono;

using option::Option;
using option::Descriptor;
using option::Parser;
using option::Stats;
using option::ArgStatus;


struct Arg: public option::Arg {
  static void printError(const char* msg1, const option::Option& opt, const char* msg2) {
    fprintf(stderr, "%s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }
  static option::ArgStatus Unknown(const option::Option& option, bool msg) {
    if (msg) printError("Unknown option '", option, "'\n");
    return option::ARG_ILLEGAL;
  }
  static option::ArgStatus NonEmpty(const option::Option& option, bool msg) {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  }
  static option::ArgStatus Numeric(const option::Option& option, bool msg) {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;
    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
};


enum  optionIndex {
  UNKNOWN,
  LIST,
  KMER,
  FETCH,
  WORD,
  HHL,
  QUERY,
  OUTPUT,
  LOGO,
  HELP
};



const option::Descriptor usage[] = {
  {UNKNOWN, 0,"" , "" , Arg::Unknown,0},
  {LIST, 0, "l" , "list" ,Arg::NonEmpty,
    "  --list, -l <file> "
      "\tUse the content of the given (raw formatted) <file> to load genomes.\v"
      "Example:"
      "\v  --list my_genomes.txt" },
  {KMER,  0, "K" , "kmer"  ,Arg::Numeric,
    "  --kmer, -K <int> "
      "\tSet the value of paramter K to the given value.\v"
      "Example:"
      "\v  --kmer 31 or -K 31" },
  {FETCH,  0, "F" , "Fetch"  ,Arg::Numeric,
    "  --Fetch, -F <int> "
      "\tSet the value of paramter F to the given value.\v"
      "Example:"
      "\v  --Fetch 16 or -F 16" },
  {WORD,  0, "W" , "Word"  ,Arg::Numeric,
    "  --Word, -W <int> "
      "\tSet the value of paramter W to the given value.\v"
      "Example:"
      "\v  --Word 10 or -W 10" },
  {HHL,  0, "H" , "HHL"  ,Arg::Numeric,
    "  --HHL, -H <int> "
      "\tSet the value of paramter H to the given value.\v"
      "Example:"
      "\v  --HHL 4 or -H 4" },
  {QUERY, 0, "Q", "query"    , Arg::NonEmpty,
    "  --query, -Q <filename> "
      "\tFor each sequence in the <sequence_files> search the sequence in the index and print the genomes"
      " with this sequence.\v"
      "The query file can be either a fasta formatted file (each sequence being a query) or a one line "
      "Examples:"
      "\v --query 'sequence_files.txt'" },
  {OUTPUT, 0, "o", "output", Arg::NonEmpty,
    "  --output, -o <filename> "
      "\tDump the current index in text format to the given file."
  },
  {LOGO, 0, "",  "logo", Arg::None,
    "  --logo "
      "\tPrint ASCII art logo, then exit."
  },
  {UNKNOWN, 0,"" , ""    , Arg::Unknown, "\n* Other usage:"},
  {HELP,  0, "h" , "help"  ,Arg::None,
    "  --help, -h  "
      "\tPrint usage and exit." },
  {0,0,0,0,0,0}
};


option::Option *options = NULL, *buffer = NULL;
void deleteOptsArrays() {
  if (options) {
    delete [] options;
  }
  if (buffer) {
    delete [] buffer;
  }
}



int main(int argc, char * argv[]){
  int F=16,K=31,W=10,H=4;
  string list_file = "";
  string query_file = "";
  string out_file = "";
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats stats(true, usage, argc, argv);
  options = new option::Option[stats.options_max];
  buffer = new option::Option[stats.buffer_max];
  atexit(deleteOptsArrays);
  option::Parser parse(usage, argc, argv, options, buffer);

  if (parse.error()) {
    cout << "Bad usage!!!" << endl;
    return EXIT_FAILURE;
  }


  /**********************************/
  /* Check Help             options */
  /**********************************/
  if (options[HELP] || argc == 0) {
    option::printUsage(clog, usage);
    return EXIT_SUCCESS;
  }

  /***********************************/
  /* Set the k-mer length and Other  */
  /***********************************/
  if (options[KMER]) {
    K = atoi(options[KMER].last()->arg);
    cout << "K,F,H,W = " <<K <<","<< F <<"," << H <<","<< W << endl;

  }
  if (options[F]) {
    F = atoi(options[FETCH].last()->arg);
    cout << "K,F,H,W = " <<K <<","<< F <<"," << H <<","<< W << endl;


  }
  if (options[H]) {
    H = atoi(options[HHL].last()->arg);
    cout << "K,F,H,W = " <<K <<","<< F <<"," << H <<","<< W << endl;
  }
  if (options[W]) {
    W = atoi(options[WORD].last()->arg);
    cout << "K,F,H,W = " <<K <<","<< F <<"," << H <<","<< W << endl;
  }

  /************************************/
  /* Complain about unknown arguments */
  /************************************/
  for (int i = 0; i < parse.nonOptionsCount(); ++i) {
    cout << "Non-option argument #" << i << " is " << parse.nonOption(i)<<endl;
    cout << "Ignoring unknown argument '" << parse.nonOption(i) << "'" << endl;
  }
  if (parse.nonOptionsCount()) {
    cout << "Bad usage!!!" << endl;
    return EXIT_FAILURE;
  }

  cout << "K,F,H,W = " <<K <<","<< F <<"," << H <<","<< W << endl;
  //    Index(uint32_t lF, uint32_t K, uint32_t W, uint32_t H);
  Index monidex(F,K,W,H);
  time_point<system_clock> start, endindex,end;
  start = std::chrono::system_clock::now();

  /*****************************************/
  /* Add the genomes given in config files */
  /*****************************************/
  if (options[LIST]) {
    list_file = options[LIST].last()->arg;       
    ifstream ifs(list_file);
    if (!ifs) {
      cout << "Unable to open the file '" << list_file << "'" << endl;
    }
    DEBUG_MSG("Opening file : '"<<list_file<<"'");
    monidex.insert_file_of_file_whole(list_file);
    DEBUG_MSG("File added");
  }
  endindex = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = endindex - start;
  std::cout << "Indexing lasted " << elapsed_seconds.count() << "s\n";

  /*****************************************/
  /* Add the query file and do the request */
  /*****************************************/
  if (options[QUERY]) {
    query_file = options[QUERY].last()->arg;       
    ifstream ifs(query_file);
    if (!ifs) {
      cout << "Unable to open the file '" << query_file << "'" << endl;
    }
    DEBUG_MSG("Opening file...");
    monidex.query_file_of_file_whole(query_file);
    DEBUG_MSG("Query done.");
  }
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - endindex;
  cout << "Query lasted " << elapsed_seconds.count() << "s\n";
  elapsed_seconds = end - start;
  cout<<"whole run tool took " << elapsed_seconds.count() << endl;


    if (options[OUTPUT]) {
      out_file = options[OUTPUT] ? (options[OUTPUT].last()->arg) : "nihmOutput";
      DEBUG_MSG("Output file name = " << out_file);
      monidex.toFile(out_file);
    }

  /**********************************************/
  /* Display the ASCII art logo of the program. */
  /**********************************************/
  if (options[LOGO]) {
    string logo_name="../resources/nihm.ascii";    
    ifstream logo;
    string line;
    logo.open(logo_name);
    if (logo.is_open()){
      while ( getline (logo,line) ){
        cout << line << '\n';
      }
      logo.close();
    }  
    else cout << "Unable to open file :'"<<logo_name<<"'"<<endl;
    return EXIT_SUCCESS;
  }


  return EXIT_SUCCESS;
}
