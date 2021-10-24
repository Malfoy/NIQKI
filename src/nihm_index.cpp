#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "nihm_index.h"
#include "common.h"



using namespace std;
const int bufferSize = 10000;



Index::Index(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iH=4, const string ifilename="nihmOutput.gz", double min_fract=0.1) {
  filename=ifilename;
  pretty_printing=false;
  lF=ilF;
  K=iK;
  W=iW;
  H=iH;
  F=1<<lF;
  M=W-H;
  min_score=min_fract*F;
  fingerprint_range=1<<W;
  mask_M=(1<<M)-1;
  maximal_remainder=(1<<H)-1;
  genome_numbers=0;
  Buckets=new vector<gid>[fingerprint_range*F];
  offsetUpdatekmer=1;
  offsetUpdatekmer<<=2*K;
  for(uint32_t i=0; i<mutex_number;++i) {
    omp_init_lock(&lock[i]);
  }
  if(pretty_printing){
      outfile=new zstr::ofstream(filename);
  }else{
      outfile=new zstr::ofstream(filename,ios::binary);
  }
}



//TODO UPDATE
void Index::dump_index_disk(const string& filestr)const{
  zstr::ofstream out(filestr);
	// VARIOUS INTEGERS
	out.write(reinterpret_cast<const char*>(&lF), sizeof(lF));
  out.write(reinterpret_cast<const char*>(&K), sizeof(K));
  out.write(reinterpret_cast<const char*>(&H), sizeof(H));
  out.write(reinterpret_cast<const char*>(&W), sizeof(W));
  out.write(reinterpret_cast<const char*>(&genome_numbers), sizeof(genome_numbers));
  for(uint i(0);i<fingerprint_range*F;++i){
    uint32_t size(Buckets[i].size());
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    out.write(reinterpret_cast<const char*>(&Buckets[i]), size*sizeof(gid));
  }
  for(uint i(0);i<filenames.size();++i){
    out.write((filenames[i]+'\n').c_str(),filenames[i].size()+1);
  }
}



//TODO UPDATE
Index::Index(const string& filestr) {
  zstr::ifstream in(filestr);
  in.read(reinterpret_cast< char*>(&lF), sizeof(lF));
  in.read(reinterpret_cast< char*>(&K), sizeof(K));
  in.read(reinterpret_cast< char*>(&H), sizeof(H));
  in.read(reinterpret_cast< char*>(&W), sizeof(W));
  in.read(reinterpret_cast< char*>(&genome_numbers), sizeof(genome_numbers));
  F=1<<lF;
  M=W-H;
  fingerprint_range=1<<W;
  mask_M=(1<<M)-1;
  maximal_remainder=(1<<H)-1;
  Buckets=new vector<gid>[fingerprint_range*F];
  offsetUpdatekmer=1;
  offsetUpdatekmer<<=2*K;
  for(uint32_t i=0; i<mutex_number;++i) {
    omp_init_lock(&lock[i]);
  }
  for(uint i(0);i<fingerprint_range*F;++i){
    uint32_t size;
    in.read(reinterpret_cast< char*>(&size), sizeof(size));
    if(size!=0){
      Buckets[i].resize(size);
      in.read(reinterpret_cast< char*>(&Buckets[i]), size*sizeof(gid));
    }
  }
}



Index::~Index() {
  delete[] Buckets;
  outfile->close();
}


//UTILS
uint64_t Index::nuc2int(char c)const {
  switch(c){
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 0;
  }
  exit(0);
  return 0;
}


string Index::kmer2str(uint64_t num)const {
  string res;
  uint64_t anc(1);
  anc<<=(2*(K-1));
  for(uint i(0);i<K;++i){
    uint nuc=num/anc;
    num=num%anc;
    if(nuc==3){
      res+="T";
    }
    if(nuc==2){
      res+="G";
    }
    if(nuc==1){
      res+="C";
    }
    if(nuc==0){
      res+="A";
    }
    if (nuc>=4){
      cout<<nuc<<endl;
      cout<<"WTF"<<endl;
    }
    anc>>=2;
  }
  return res;
}


//UTILS
uint64_t Index::asm_log2(const uint64_t x) const {
  uint64_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
      );
  return y;
}


//UTILS
uint64_t Index::nuc2intrc(char c)const {
  switch(c) {
    case 'A': return 3;
    case 'C': return 2;
    case 'G': return 1;
    default: return 0;
  }
  //~ cout<<"Unknow nucleotide: "<<c<<"!"<<endl;
  exit(0);
  return 0;
}



void Index::update_kmer(kmer& min, char nuc)const {
  min<<=2;
  min+=nuc2int(nuc);
  min%=offsetUpdatekmer;
}



void Index::update_kmer_RC(kmer& min, char nuc)const {
  min>>=2;
  min+=(nuc2intrc(nuc)<<(2*K-2));
}



kmer Index::rcb(kmer min)const {
  kmer res(0);
  kmer offset(1);
  offset<<=(2*K-2);
  for(uint i(0); i<K;++i){
    res+=(3-(min%4))*offset;
    min>>=2;
    offset>>=2;
  }
  return res;
}


//UTILS
kmer Index::str2numstrand(const string& str)const {
  uint64_t res(0);
  for(uint i(0);i<str.size();i++) {
    res<<=2;
    switch (str[i]){
      case 'A':res+=0;break;
      case 'C':res+=1;break;
      case 'G':res+=2;break;
      case 'T':res+=3;break;

      case 'a':res+=0;break;
      case 'c':res+=1;break;
      case 'g':res+=2;break;
      case 't':res+=3;break;
      default: return 0 ;break;
    }
  }
  return (res);
}



int32_t Index::get_fingerprint(uint64_t hashed)const {
	//~ print_bin(hashed);
  int32_t result;
  result = hashed&mask_M;//we keep the last bits for the minhash part
  uint32_t ll=asm_log2(hashed);//we compute the log of the hash
   //~ cout<<ll<<endl;
  uint32_t size_zero_trail(64-ll-1);
  //~ cout<<size_zero_trail<<endl;
  int remaining_nonzero=maximal_remainder-size_zero_trail;
  //~ cout<<remaining_nonzero<<endl;
  remaining_nonzero=max(0,remaining_nonzero);
  
   //~ cout<<remaining_nonzero<<endl;
   //~ cin.get();
  // if the log is too high can cannont store it on H bit here we sature
  result+=remaining_nonzero<<M;// we concatenant the hyperloglog part with the minhash part
  //~ print_bin(result,W);
  //~ cin.get();
  return result;
}



uint64_t Index::revhash64 ( uint64_t x ) const {
  x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
  x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
  x = ( ( x >> 32 ) ^ x );
  return x;
}



uint64_t unrevhash64(uint64_t x) {
  x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
  x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
  x = ((x >> 32) ^ x);
  return x;
}



void Index::compute_sketch(const string& reference, vector<int32_t>& sketch) const {
  if(sketch.size()!=F) {
    sketch.resize(F,-1);
  }
  kmer S_kmer(Index::str2numstrand(reference.substr(0,K-1)));//get the first kmer (k-1 bases)
  kmer RC_kmer(Index::rcb(S_kmer));//The reverse complement
  for(uint i(0);i+K<reference.size();++i) {// all kmer in the genome
    Index::update_kmer(S_kmer,reference[i+K-1]);
    Index::update_kmer_RC(RC_kmer,reference[i+K-1]);
    kmer canon(min(S_kmer,RC_kmer));//Kmer min, the one selected
    uint64_t hashed=revhash64(canon);
    uint64_t bucket_id(unrevhash64(canon)>>(64-lF));//Which Bucket 
    int32_t fp=get_fingerprint(hashed);
    //MAXHASH
	//~ if(sketch[bucket_id] < fp or sketch[bucket_id]<0) {
	//MINHASH
    if(sketch[bucket_id] > fp or sketch[bucket_id]==-1) {
      sketch[bucket_id]=fp;
    }
  }
}



//HERE we only select the minimal hashes without computing the HMH fingerprint
//~ void Index::compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch)const{
  //~ if(sketch.size()!=F) {
    //~ sketch.resize(F,-1);
  //~ }
  //~ kmer S_kmer(Index::str2numstrand(reference.substr(0,K-1)));
  //~ kmer RC_kmer(Index::rcb(S_kmer));
  //~ for(uint i(0);i+K<reference.size();++i) {
    //~ Index::update_kmer(S_kmer,reference[i+K-1]);
    //~ Index::update_kmer_RC(RC_kmer,reference[i+K-1]);
    //~ kmer canon(min(S_kmer,RC_kmer));
    //~ uint64_t hashed=revhash64(canon);
    //~ uint32_t bucket_id(unrevhash64(canon)>>(64-lF));
    //~ if((uint64_t)sketch[bucket_id]>hashed || sketch[bucket_id] == (uint64_t)-1) {
      //~ sketch[bucket_id]=hashed;
    //~ }
  //~ }
//~ }



void Index::insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id) {
  for(int i(0);i<F;++i) {
    if(sketch[i]<fingerprint_range and sketch[i]>=(int)0) {
      omp_set_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
      Buckets[sketch[i]+i*fingerprint_range].push_back(genome_id);
      omp_unset_lock(&lock[(sketch[i]+i*fingerprint_range)%mutex_number]);
    }
  }
}



void Index::insert_sequence(const string& str,uint32_t genome_id) {
  vector<int32_t> sketch;
  compute_sketch(str,sketch);
  insert_sketch(sketch,genome_id);
}



//all the lines from the file is considered as a separate entry with a different identifier
void Index::insert_file_lines(const string& filestr) {
  char type=get_data_type(filestr);
  zstr::ifstream in(filestr);
  #pragma omp parallel
  {
    string ref,header;
    uint32_t id;
    while(not in.eof()) {
      #pragma omp critical (in)
      {
        Biogetline(&in,ref,type,header);
      }
      if(ref.size()>K) {
        #pragma omp critical (genome_numbers)
      {
        id=genome_numbers;
        genome_numbers++;
        filenames.push_back(header);
      }
        insert_sequence(ref,id);

      }
      ref.clear();
    }
  }
}



void Index::query_file_lines(const string& filestr)const {
  char type=get_data_type(filestr);
  zstr::ifstream in(filestr);
   #pragma omp parallel
   {
    string ref,head;
    while(not in.eof()) {
      #pragma omp critical (in)
      {
        Biogetline(&in,ref,type);
      }
      if(ref.size()>K) {
        auto out(query_sequence(ref));
        ref.clear();
        output_query(out,filestr);
      }
    }
   }
}



void Index::merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2) const {
  for(uint i(0);i<sketch1.size();++i) {
    sketch1[i]=min(sketch1[i],sketch2[i]);
  }
}



void Index::insert_file_whole(const string& filestr,uint32_t identifier) {
  char type=get_data_type(filestr);
  zstr::ifstream in(filestr);
  string ref;
  vector<uint64_t> kmer_sketch;
  vector<int32_t> sketch(F,-1);
  while(not in.eof()) {
    Biogetline(&in,ref,type);
    if(ref.size()>=K) {
        compute_sketch(ref,sketch);
    }
    ref.clear();
  }   
  insert_sketch(sketch,identifier);
}



//HERE all the files of the fof are inserted as a separate entry in the index
void Index::insert_file_of_file_whole(const string& filestr) {
  DEBUG_MSG("Opening file : '"<<filestr<<"'");
  ifstream in(filestr);
  if (!in) {
    cout << "Unable to open the file '" << filestr << "'" << endl;
    exit(0);
  }
//~ #pragma omp parallel
{
    string ref;
    uint32_t id;
    while(not in.eof()) {
      bool go=false;
      #pragma omp critical
      {
        getline(in,ref);
        DEBUG_MSG("Getline from file :'"<<filestr<<"' = '"<<ref<<"'");

        if(ref.size()>2){
          if(exists_test(ref)) {
            id=genome_numbers;
            DEBUG_MSG("Genome numbers: "<<genome_numbers);
            genome_numbers++;
            filenames.push_back(ref);
            go=true;
          }
        }
      }
      if(go) {
        DEBUG_MSG("Adding file :'"<<ref<<"'");
        insert_file_whole(ref,id);
        DEBUG_MSG("File: '"<<ref<<"' added");
      }
    }
  }
}



//HERE all the kmer of the file are put in a single sketch and Queried
void Index::query_file_whole(const string& filestr) {
  char type=get_data_type(filestr);
  zstr::ifstream in(filestr);
  string ref;
  vector<int32_t> sketch(F,-1);
  while(not in.eof()){
    Biogetline(&in,ref,type);
    if(ref.size()>K){
      compute_sketch(ref,sketch);
    }

  }   

  auto out(query_sketch(sketch));
  output_query(out,filestr);
}
//~ atomic<uint64_t> sum41(0),denom(1);



void Index::query_file_of_file_whole(const string& filestr) {
  zstr::ifstream in(filestr);
//~ #pragma omp parallel
{
  string ref;
  while(not in.eof()){
    
  #pragma omp critical (input)
    {
      getline(in,ref);
    }
    if(exists_test(ref)){
      query_file_whole(ref);
    }
    ref.clear();
  }
}
//~ cout<<sum41/denom<<endl;
}




void Index::output_query(const query_output& toprint,const string& queryname)const{
  if(pretty_printing){
    #pragma omp critical (outputfile)
    {
		
      *outfile<<queryname<<"\n";
      for(uint i(0);i<toprint.size();++i){
        *outfile<<filenames[toprint[i].second]<<" "<<toprint[i].first<<'\n';
        //~ if(toprint[i].first!=F){
			//~ sum41+=toprint[i].first;
			//~ denom++;
		//~ }
      }
      
    }
  }else{
    #pragma omp critical (outputfile)
    {
      *outfile<<queryname<<"\n";
      uint32_t size(toprint.size());
      outfile->write(reinterpret_cast<const char*>(&(size)), sizeof(size));
      for(uint i(0);i<toprint.size();++i){
		 outfile->write(reinterpret_cast<const char*>(&(toprint[i].second)), sizeof(toprint[i].second));
		 outfile->write(reinterpret_cast<const char*>(&(toprint[i].first)), sizeof(toprint[i].first));
      }
    }
  } 
}



query_output Index::query_sketch(const vector<int32_t>& sketch)const {
    query_output result;
     if(lF<=7){
        uint8_t counts[genome_numbers]={0};
        for(int i(0);i<F;++i){
            if(sketch[i]<(int32_t)fingerprint_range and sketch[i]>=0){
                for(uint j(0);j<Buckets[sketch[i]+i*fingerprint_range].size();++j){
                  counts[Buckets[sketch[i]+i*fingerprint_range][j]]++;
                }
            }
        }
        for(uint32_t i(0);i<genome_numbers;++i){
            if((uint32_t)counts[i]>=min_score){
                result.push_back({counts[i],i});
            }
        }
      }else if(lF<=15){
        uint16_t counts[genome_numbers]={0};
        for(int i(0);i<F;++i){
            if(sketch[i]<(int32_t)fingerprint_range and sketch[i]>=0){
                for(uint j(0);j<Buckets[sketch[i]+i*fingerprint_range].size();++j){
                  counts[Buckets[sketch[i]+i*fingerprint_range][j]]++;
                }
            }
        }
        for(uint32_t i(0);i<genome_numbers;++i){
            if((uint32_t)counts[i]>=min_score){
                result.push_back({counts[i],i});
            }
        }
    }else{
        uint32_t counts[genome_numbers]={0};
        for(int i(0);i<F;++i){
            if(sketch[i]<(int32_t)fingerprint_range and sketch[i]>=0){
                for(uint j(0);j<Buckets[sketch[i]+i*fingerprint_range].size();++j){
                    counts[Buckets[sketch[i]+i*fingerprint_range][j]]++;
                }
            }
        }
        for(uint32_t i(0);i<genome_numbers;++i){
            if((uint32_t)counts[i]>=min_score){
                result.push_back({counts[i],i});
            }
        }
    }
    
    sort(result.begin(),result.end(),greater<pair<uint32_t,uint32_t>>());
    return result;
}



query_output Index::query_sequence(const string& str)const {
    vector<int32_t> sketch;
    compute_sketch(str,sketch);
    return query_sketch(sketch);
}

/*******************************************************************************************/
/******************************* Code Duplication for matrix *******************************/
/*******************************************************************************************/

void Index::query_file_whole_matrix(const string& filestr) {
  char type=get_data_type(filestr);
  zstr::ifstream in(filestr);
  string ref;
  vector<int32_t> sketch(F,-1);
  while(not in.eof()){
    Biogetline(&in,ref,type);
    if(ref.size()>K){
      compute_sketch(ref,sketch);
    }
  }

  auto out(query_sketch(sketch));
  DEBUG_MSG("Sending out in :'"<<filestr<<"'");
  output_matrix(out,filestr);
}

void Index::query_file_of_file_whole_matrix(const string& filestr) {
  //TODO LE HEADER
  *outfile<<"##Names"<<"\t";
  for(uint i(0);i<filenames.size();++i){
    *outfile<<filenames[i]<<"\t";
  }
  *outfile<<"\n";
  zstr::ifstream in(filestr);
#pragma omp parallel
{
  string ref;
  while(not in.eof()){

    #pragma omp critical (input)
    {
      getline(in,ref);
    }
    if(exists_test(ref)){
      query_file_whole_matrix(ref);
      *outfile<<endl;
    }
    ref.clear();
  }
}
}

void Index::output_matrix(const query_output& toprint,const string& queryname)const{
#pragma omp critical (outputfile)
  {
    double buffer[genome_numbers]={0};
    for(uint i(0);i<toprint.size();++i){
      buffer[toprint[i].second]=((double)toprint[i].first)/F;
    }
    *outfile<<queryname<<"\t";
    for(uint i(0);i<genome_numbers;++i){
      *outfile<<buffer[i]<<"\t";
    }
  }
}
/***********************************************************************************************/
/******************************* End Code Duplication for matrix *******************************/
/***********************************************************************************************/

atomic<uint32_t> genomes_downloaded(0);
atomic<uint64_t> bases_downloaded(0);



string get_output_exec_cmd(const char* cmd) {
  array<char, 1024*1024> buffer;
  string result;
  unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {

    result += buffer.data();
  }
  return result;
}



string get_name_ncbi(const string& str){
  uint lastposition(0);
  for(uint i(0);i<str.size()-3;++i){
    if(str[i]=='/'){
      lastposition=i;
    }
  }
  lastposition++;
  return str.substr(lastposition,str.size()-lastposition);
}



bool Index::Download_NCBI(const string& str, vector<int32_t>& sketch){
  string cmd("wget -qO - "+str+"/"+get_name_ncbi(str)+"_genomic.fna.gz | gzip -d -c -");
  array<char, 1024*1024> buffer;
  string result;
  string sequence;
  bool something_to_eat(false);
  unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
  if (!pipe) {
    return false;
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result= buffer.data();
    if(result[0]!='>'){
      if(result[result.size()-1]=='\n'){
        sequence+=result.substr(0,result.size()-1);
      }else{
        sequence+=result;
      }
    }else{
      if(sequence.size()>K){
        something_to_eat=true;
        compute_sketch(sequence,sketch);
        bases_downloaded+=sequence.size();
        sequence.clear();
      }
    }
    result.clear();
  }
  if(something_to_eat){
    uint gtemp=++genomes_downloaded;
    if(gtemp%1000==0){
       cout<<"#genomes downloaded: "<<intToString(genomes_downloaded)<<" #bases downloaded: "<<intToString(bases_downloaded)<<endl;
    }
  }
  return something_to_eat;
}



//Hacky function to download genome from ncbi
void Index::Download_NCBI_fof(const string& fofncbi){
  zstr::ifstream in(fofncbi);
#pragma omp parallel
{
  string ref;
  uint32_t id;
  while(not in.eof()){
    
    #pragma omp critical (in)
    {
      getline(in,ref);
    }
    if(ref.size()>5){
      vector<int32_t> sketch(F,-1);
      if(Download_NCBI(ref,sketch)){
        #pragma omp critical (genomenumber)
        {
          id=genome_numbers;
          genome_numbers++;
          filenames.push_back(ref);
        }
        insert_sketch(sketch,id);
      }
    }
    ref.clear();
  }
}
  cout<<"#genomes downloaded: "<<intToString(genomes_downloaded)<<" #bases downloaded: "<<intToString(bases_downloaded)<<endl;
}


//pretty printing
string Index::intToString(uint64_t n) {
  if (n < 1000) {
    return to_string(n);
  }
  string end(to_string(n % 1000));
  if (end.size() == 3) {
    return intToString(n / 1000) + "," + end;
  }
  if (end.size() == 2) {
    return intToString(n / 1000) + ",0" + end;
  }
  return intToString(n / 1000) + ",00" + end;
}



void Index::Biogetline(zstr::ifstream* in,string& result,char type,string& header)const {
  string discard;
  result.clear();
  switch(type){
    case 'Q':
      getline(*in,header);
      getline(*in,result);
      getline(*in,discard);
      getline(*in,discard);
      break;
    case 'A':
      getline(*in,header);
      char c=in->peek();
      while(c!='>' and c!=EOF){
        getline(*in,discard);
        result+=discard;
        c=in->peek();
      }
      break;
  }
  if(result.size()< K){
    result.clear();
    header.clear();
  }
}



void Index::Biogetline(zstr::ifstream* in,string& result,char type)const {
  string discard;
  result.clear();
  switch(type){
    case 'Q':
      getline(*in,discard);
      getline(*in,result);
      getline(*in,discard);
      getline(*in,discard);
      break;
    case 'A':
      getline(*in,discard);
      char c=in->peek();
      while(c!='>' and c!=EOF){
        getline(*in,discard);
        result+=discard;
        c=in->peek();
      }
      break;
  }
  if(result.size()< K){
    result.clear();
  }
}


char Index::get_data_type(const string& filename)const{
  if(filename.find(".fq")!=string::npos){
    return 'Q';
  }
  if(filename.find(".fastq")!=string::npos){
    return 'Q';
  }
  return 'A';
}
