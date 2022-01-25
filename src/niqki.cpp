#include "niqki_index.h"
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
#include <iomanip> // std::setw
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
  LISTLINES,
  QUERYLINES,
  KMER,
  FETCH,
  WORD,
  HHL,
  MIN,
  QUERY,
  OUTPUT,
  DUMP,
  LOAD,
  MATRIX,
  LOGO,
  DOWNLAD,
  PRETTY,
  GENOME_SIZE,
  HELP
};



const option::Descriptor usage[] = {
  {UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Input***"},
  {LIST, 0, "I" , "index" ,Arg::NonEmpty,
    "  --index, -I <filename> "
      "\tInput file of files to Index.\v"
     },
  
   {QUERY, 0, "Q", "query"    , Arg::NonEmpty,
    "  --query, -Q <filename> "
      "\tInput file of file to Query.\v"
      },
   {LISTLINES, 0, "i" , "indexlines" ,Arg::NonEmpty,
    "  --indexlines, -i <filename> "
      "\tQuery fa/fq file where each line is a separate entry to Index\v"
     },
     {QUERYLINES, 0, "l", "querylines"    , Arg::NonEmpty,
    "  --querylines, -q <filename> "
      "\tInput fa/fq where each line is a separate entry to Query\v"
      },
  {UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Main parameters***"},
  {KMER,  0, "K" , "kmer"  ,Arg::Numeric,
    "  --kmer, -K <int> "
      "\tKmer size (31).\v"
      },
  {FETCH,  0, "S" , "sketch"  ,Arg::Numeric,
    "  --sketch, -S <int> "
      "\tSet sketch size to 2^S (15).\v"
      },
 
 {UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Output***"},
  {OUTPUT, 0, "O", "output", Arg::NonEmpty,
    "  --output, -O <filename> "
      "\tOutput file (niqkiOutput.gz)"
  },

     {MIN,  0, "J" , "minjac"  ,Arg::NonEmpty,
    "  --minjac, -J <int> "
      "\tMinimal jaccard Index to report (0.1).\v"
      },
   {PRETTY, 0, "P", "pretty", Arg::None,
    "  --pretty, -P "
      "\t Print a human-readable outfile. By default the outfile is in binary."
  },
  {MATRIX, 0, "M", "matrix", Arg::NonEmpty,
    "  --matrix, -M <filename> "
      "\tOutput the matrix distance to the given file."
  },
  {UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Advanced parameters*** (You know what you are doing)"},
   {WORD,  0, "W" , "word"  ,Arg::Numeric,
    "  --word, -W <int> "
      "\tFingerprint size (12). Modify with caution, larger fingerprints enable queries with less false positive but increase EXPONENTIALY the overhead as the index count S*2^W cells. \v"
     },
      {GENOME_SIZE,  0, "G" , "Genomes_sizes"  ,Arg::Numeric,
    "  --Genomes_sizes, -G <int> "
      "\tRought expectation of the genome sizes. \v"
      },
  {HHL,  0, "H" , "HHL"  ,Arg::Numeric,
    "  --HHL, -H <int> "
      "\tSize of the hyperloglog section (4).  Modify with caution and prefer to use -G.\v"
      },
{UNKNOWN, 0,"" , "" , Arg::Unknown,"\n***Index files***"},
  {DUMP, 0, "D", "dump", Arg::NonEmpty,
    "  --dump, -D <filename> "
      "\tDump the current index to the given file."
  },
  {LOAD, 0, "L", "load", Arg::NonEmpty,
    "  --load, -L <filename> "
      "\tLoad an index to the given file."
  },
    {UNKNOWN, 0,"" , ""    , Arg::Unknown, "\n***Other***"},
    {DOWNLAD, 0, "Iddl" , "indexdownload" ,Arg::NonEmpty,
    "  --indexdownload, -Iddl <filename> "
      "\tGet a list of NCBI accesion to download and to put it in the index (experimental). This this post to get such a list: https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/#allcomplete \v"
     },
  {LOGO, 0, "",  "logo", Arg::None,
    "  --logo "
      "\tPrint ASCII art logo, then exit."
  },

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



static int old_wd;
static void changeDirFromFilename(const char* fname) {
  DEBUG_MSG("CWD is " << getcwd(NULL, 0));
  old_wd = open(".", O_CLOEXEC);
  char fname_copy[PATH_MAX];
  strncpy(fname_copy, fname, PATH_MAX);
  fname_copy[PATH_MAX - 1] = '\0';
  char *dname = dirname(fname_copy);
  DEBUG_MSG("Changing to directory " << dname);
  errno = 0;
  if (chdir(dname)) {
    cout << "Error: " << strerror(errno) << endl;
  }
  DEBUG_MSG("Now CWD is " << getcwd(NULL, 0));
}



static void restoreDir() {
  errno = 0;
  if (fchdir(old_wd)) {
    cout << "Error: " << strerror(errno) << endl;
  }
  DEBUG_MSG("Restore working directory to " << getcwd(NULL, 0));
}



int main(int argc, char * argv[]){
  int F=0,K=0,W=0,H=0;
  bool pretty_printing=true;
  double min_fract;
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
    option::printUsage(clog, usage,160);
    return EXIT_SUCCESS;
  }

  /***********************************/
  /* Set the k-mer length and Other  */
  /***********************************/
  K = options[KMER] ? atoi(options[KMER].last()->arg) : 31;
  DEBUG_MSG("K = " << K);
  F = options[FETCH] ? atoi(options[FETCH].last()->arg) : 15;
  DEBUG_MSG("F = " << F);
  H = options[HHL] ? atoi(options[HHL].last()->arg) : 4;
  DEBUG_MSG("H = " << H);
  W = options[WORD] ? atoi(options[WORD].last()->arg) : 12;
  DEBUG_MSG("W = " << W);
  min_fract = options[MIN] ? atof(options[MIN].last()->arg) : 0;
  DEBUG_MSG("min_fract = " << min_fract);
  uint genomes_sizes = options[GENOME_SIZE] ? atoi(options[GENOME_SIZE].last()->arg) : 0;

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
  if (options[OUTPUT]) {
    out_file = options[OUTPUT] ? (options[OUTPUT].last()->arg) : "niqkiOutput.gz";
    DEBUG_MSG("Output file name = " << out_file);
  } else {
    out_file="niqkiOutput.gz";
  }
  cout << "+-------------------------------------------------------------------+" << endl;
  cout << "|                            Informations                           |" << endl;
  cout << "+-----------------------------------+-------------------------------+" << endl;
  Index* monindex;
   if (options[PRETTY]){
    pretty_printing=true;
  } 
  if(options[LOAD]){
    string indexfile = options[LOAD].last()->arg;  
    monindex = new Index(indexfile,pretty_printing,out_file);
  }else{
    monindex=new Index(F,K,W,H,out_file,min_fract);
    monindex->pretty_printing=pretty_printing;
  }
   if(genomes_sizes!=0){
    monindex->select_best_H(genomes_sizes);
  }
 
  time_point<system_clock> start, endindex,end;
  start = std::chrono::system_clock::now();

  /*****************************************/
  /* Add the genomes given in config files */
  /*****************************************/
  if (options[LIST]) {
    //option::Option* opt = options[LIST];
    list_file = options[LIST].last()->arg;       
    ifstream ifs(list_file);
    if (!ifs) {
      cout << "Unable to open the file '" << list_file << "'" << endl;
    }
    changeDirFromFilename(list_file.c_str());
    DEBUG_MSG("Opening file : '"<<list_file<<"'");
    monindex->insert_file_of_file_whole(list_file.substr(list_file.find_last_of("/\\") + 1));
    DEBUG_MSG("File added");
    restoreDir();

  }

  if (options[LISTLINES]) {
    //option::Option* opt = options[LIST];
    list_file = options[LISTLINES].last()->arg;       
    ifstream ifs(list_file);
    if (!ifs) {
      cout << "Unable to open the file '" << list_file << "'" << endl;
    }
    changeDirFromFilename(list_file.c_str());
    DEBUG_MSG("Opening file : '"<<list_file<<"'");
    monindex->insert_file_lines(list_file.substr(list_file.find_last_of("/\\") + 1));
    DEBUG_MSG("File added");
    restoreDir();
  }


  if (options[DOWNLAD]) {
    //option::Option* opt = options[LIST];
    list_file = options[DOWNLAD].last()->arg;       
    ifstream ifs(list_file);
    if (!ifs) {
      cout << "Unable to open the file '" << list_file << "'" << endl;
    }
    changeDirFromFilename(list_file.c_str());
    DEBUG_MSG("Opening file : '"<<list_file<<"'");
    monindex->Download_NCBI_fof(list_file.substr(list_file.find_last_of("/\\") + 1));
    DEBUG_MSG("File added");
    restoreDir();
  }

  if(options[DUMP]){
    string indexfile = options[DUMP].last()->arg;  
    monindex->dump_index_disk(indexfile);
  }


  endindex = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = endindex - start;
  cout << "| Indexing lasted (s)               |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;


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
    monindex->query_file_of_file_whole(query_file);
    DEBUG_MSG("Query done.");
  }

    if (options[QUERYLINES]) {
    query_file = options[QUERYLINES].last()->arg;       
    ifstream ifs(query_file);
    if (!ifs) {
      cout << "Unable to open the file '" << query_file << "'" << endl;
    }
    DEBUG_MSG("Opening file...");
    monindex->query_file_lines(query_file);
    DEBUG_MSG("Query done.");
  }
  monindex->outfile->close();


  end = std::chrono::system_clock::now();
  elapsed_seconds = end - endindex;
  cout << "| Query lasted (s)                  |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;
  elapsed_seconds = end - start;
  cout << "| Whole run lasted (s)              |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;

    if (options[MATRIX]) {
    string matrix_file = options[MATRIX].last()->arg;
    ifstream ifs(matrix_file);
    if (!ifs) {
      cout << "Unable to open the file '" << matrix_file << "'" << endl;
    }
    if (!options[LIST] and !options[LISTLINES]) {
	  start= std::chrono::system_clock::now();
      changeDirFromFilename(matrix_file.c_str());
      DEBUG_MSG("Creating the index from: '"<<matrix_file.substr(matrix_file.find_last_of("/\\") + 1)<<"'");
      monindex->insert_file_of_file_whole(matrix_file.substr(matrix_file.find_last_of("/\\") + 1));
      DEBUG_MSG("File added");
      restoreDir();
      DEBUG_MSG("Directory restored");
      endindex = std::chrono::system_clock::now();
	  std::chrono::duration<double> elapsed_seconds = endindex - start;
	  cout << "| Indexing lasted (s)               |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;
    }
    changeDirFromFilename(matrix_file.c_str());
    start= std::chrono::system_clock::now();
    monindex->query_matrix();
    end= std::chrono::system_clock::now();
    elapsed_seconds = end - start;
	cout << "| Query lasted (s)                  |" << setw(30) << setfill(' ') << elapsed_seconds.count() << " |" << endl;
    restoreDir();
  }


  /**********************************************/
  /* Display the ASCII art logo of the program. */
  /**********************************************/
  if (options[LOGO]) {
    string logo_name="../resources/niqki.ascii";    
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
  cout << "+-----------------------------------+-------------------------------+" << endl;
  cout << "| k-mer size                        |" << setw(30) << setfill(' ') << K << " |" << endl
       << "| S                                 |" << setw(30) << setfill(' ') << F << " |" << endl
       << "| Number of fingerprints            |" << setw(30) << setfill(' ') << monindex->F<< " |" << endl
       << "| W                                 |" << setw(30) << setfill(' ') << W << " |" << endl
       << "| H                                 |" << setw(30) << setfill(' ') << H << " |" << endl
       << "| Number of indexed genomes         |" << setw(30) << setfill(' ') << monindex->getNbGenomes() << " |" << endl;
  cout << "+-----------------------------------+-------------------------------+" << endl;

  return EXIT_SUCCESS;
}
