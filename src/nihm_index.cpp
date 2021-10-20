#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "nihm_index.h"
#include "common.h"


using namespace std;

Index::Index(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iH=4) {
  lF=ilF;
  K=iK;
  W=iW;
  H=iH;
  F=1<<lF;
  M=W-H;
  fingerprint_range=1<<W;
  mask_M=(1<<M)-1;
  maximal_remainder=(1<<H)-1;
  genome_numbers=0;
  mask_fingerprint=((uint64_t)1<<(64-lF))-1;
  Buckets=new vector<gid>[fingerprint_range*F];
  offsetUpdatekmer=1;
  offsetUpdatekmer<<=2*K;
  for(uint32_t i=0; i<mutex_number;++i) {
    omp_init_lock(&lock[i]);
  }
}



Index::~Index() {
  delete[] Buckets;
}



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



string Index::kmer2str(uint64_t num,uint k)const {
  string res;
  uint64_t anc(1);
  anc<<=(2*(k-1));
  for(uint i(0);i<k;++i){
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




uint64_t Index::asm_log2(const uint64_t x) const {
  uint64_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
      );
  return y;
}


uint64_t Index::nuc2intrc(char c)const {
  switch(c) {
    /*
       case 'a': return 0;
       case 'c': return 1;
       case 'g': return 2;
       case 't': return 3;
     */
    case 'A': return 3;
    case 'C': return 2;
    case 'G': return 1;
              //~ case 'T': return 0;
              //~ case 'N': return 0;
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



uint32_t Index::get_fingerprint(uint64_t hashed)const {
  // cout<<"get_fingerprint"<<endl;
  // cout<<hashed<<endl;
  // print_bin(hashed,64-lF);
  uint32_t result;
  result = hashed&mask_M;//we keep the last bits for the minhash part
  uint32_t ll=asm_log2(hashed);//we compute the log of the hash
  // cout<<ll<<endl;
  // cout<<maximal_remainder<<endl;
  uint32_t size_zero_trail(64-lF-ll-1);
  // cout<<size_zero_trail<<endl;
  int remaining_nonzero=maximal_remainder-size_zero_trail;
  remaining_nonzero=max(0,remaining_nonzero);
  // if the log is too high can cannont store it on H bit here we sature
  // cout<<ll<<endl;
  // print_bin(result,8);
  result+=remaining_nonzero<<M;// we concatenant the hyperloglog part with the minhash part
  // cout<<"get_fingerprintOUT"<<endl;
  // cout<<remaining_nonzero<<endl;
  // print_bin(result,8);
  // cin.get();
  return result;
}



uint64_t Index::revhash64 ( uint64_t x ) const {
  x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
  x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
  x = ( ( x >> 32 ) ^ x );
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
    uint32_t bucket_id(hashed>>(64-lF));//Which Bucket 
    // print_bin(hashed);
    hashed&=mask_fingerprint;
    // print_bin(hashed);
    // print_bin(sketch[bucket_id]);
    hashed=get_fingerprint(hashed);
    if(sketch[bucket_id]<hashed || sketch[bucket_id] == -1) {
      sketch[bucket_id]=hashed;
    }
  }
}



//HERE we only select the minimal hashes without computing the HMH fingerprint
void Index::compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch)const{
  // cout<<"compute_sketch_kmer"<<endl;
  if(sketch.size()!=F) {
    sketch.resize(F,-1);
  }
  kmer S_kmer(Index::str2numstrand(reference.substr(0,K-1)));
  kmer RC_kmer(Index::rcb(S_kmer));
  for(uint i(0);i+K<reference.size();++i) {
    Index::update_kmer(S_kmer,reference[i+K-1]);
    Index::update_kmer_RC(RC_kmer,reference[i+K-1]);
    kmer canon(min(S_kmer,RC_kmer));
    uint64_t hashed=revhash64(canon);
    uint32_t bucket_id(hashed>>(64-lF));
    hashed&=mask_fingerprint;
    // cout<<bucket;_id<<" "<<sketch.size()<<endl;
    if(sketch[bucket_id]<hashed || sketch[bucket_id] == -1) {
      sketch[bucket_id]=hashed;
    }
  }
  // cout<<"compute_sketch_kmer end"<<endl;
}



void Index::insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id) {
      for(uint i(0);i<F;++i) {
          if(sketch[i]<fingerprint_range and sketch[i]>=0) {
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
//TODO HANDLE FASTQ multiFASTA
void Index::insert_file_lines(const string& filestr) {
  zstr::ifstream in(filestr);
  string ref,head;
  while(not in.eof()) {
    {
      getline(in,head);
      getline(in,ref);
    }
    if(ref.size()>K) {
      insert_sequence(ref,genome_numbers);
      genome_numbers++;
    }
    ref.clear();

  }
}



void Index::query_file_lines(const string& filestr) {
  zstr::ifstream in(filestr);
  string ref,head;
  while(not in.eof()) {
    {
      getline(in,head);
      getline(in,ref);
    }
    if(ref.size()>K) {
      auto out(query_sequence(ref));
      ref.clear();
      cout<<"out: "<<endl;
      for(uint32_t i=0;i<min((uint)5,(uint)out.size());i++) {
        cout<<out[i].second<<" "<<out[i].first<<endl;
      }
    }

  }
}



void Index::merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2) {
  for(uint i(0);i<sketch1.size();++i) {
    sketch1[i]=min(sketch1[i],sketch2[i]);
  }
}



//HERE all the kmer of the file are put in a single sketch and inserted
void Index::insert_file_whole(const string& filestr) {
  zstr::ifstream in(filestr);
  string ref,head;
  vector<uint64_t> kmer_sketch;
  vector<int32_t> sketch(F,-1);
  while(not in.eof()) {
      getline(in,head);
      getline(in,ref);
    if(ref.size()>K) {
      //  cout<<"get compute_sketch_kmer"<<endl;
      compute_sketch_kmer(ref,kmer_sketch);
      // cout<<"get compute_sketch_kmer END"<<endl;
    }else{
      cout<<ref<<endl;
    }
    ref.clear();
  }   
  // cout<<"get fingerprint"<<endl;
  for(uint i(0);i<F;++i) {
    sketch[i]=get_fingerprint(kmer_sketch[i]);
  }
  // cout<<"get fingerprintOK"<<endl;
  insert_sketch(sketch,genome_numbers);
  // cout<<"Insertion done"<<endl;
  genome_numbers++;
}


void Index::insert_file_whole(const string& filestr,uint32_t identifier) {
  zstr::ifstream in(filestr);
  string ref,head;
  vector<uint64_t> kmer_sketch;
  vector<int32_t> sketch(F,-1);
  while(not in.eof()) {
  getline(in,head);
  getline(in,ref);
  if(ref.size()>K) {
      compute_sketch_kmer(ref,kmer_sketch);
  }
  ref.clear();
  }   
  // cout<<"get fingerprint"<<endl;
  for(uint i(0);i<F;++i) {
      sketch[i]=get_fingerprint(kmer_sketch[i]);
  }
  // cout<<"get fingerprintOK"<<endl;
  insert_sketch(sketch,identifier);
}



//HERE all the files of the fof are inserted as a separate entry in the index
    void Index::insert_file_of_file_whole(const string& filestr) {
        ifstream in(filestr);
        // cout<<"insert_file_of_file_whole GO"<<endl;
        //DEBUG_MSG("File name = '"<<filestr<<"'");
        #pragma omp parallel
        while(not in.eof()) {
			string ref;
            uint32_t id;
            #pragma omp critical
            {
              getline(in,ref);
              id=genome_numbers;
              genome_numbers++;
              //DEBUG_MSG("Genome numbers : "<< genome_numbers);
            }
            if(exists_test(ref)) {
              // cout<<"insert_file_whole GO"<<endl;
              insert_file_whole(ref,id);
              // cout<<"insert_file_whole END"<<endl;
            }
        }
    }



//HERE all the kmer of the file are put in a single sketch and Queried
void Index::query_file_whole(const string& filestr) {
  zstr::ifstream in(filestr);
  string ref,head;
  vector<uint64_t> kmer_sketch;
  vector<int32_t> sketch(F,-1);
  while(not in.eof()){
    {
      getline(in,head);
      getline(in,ref);
    }
    if(ref.size()>K){
      compute_sketch_kmer(ref,kmer_sketch);
    }

  }   
  for(uint i(0);i<F;++i) {
    // cout<<i<<" "<<sketch.size()<<endl;
    sketch[i]=get_fingerprint(kmer_sketch[i]);
  }
  //DEBUG_MSG("Genome Name : "<<filestr);
  // cout<<"la fingernbtoint"<<endl;
  auto out(query_sketch(sketch,10000));
  // cout<<out.size()<<endl;
}



void Index::query_file_of_file_whole(const string& filestr) {
  zstr::ifstream in(filestr);
  
  #pragma omp parallel
  while(not in.eof()){
    string ref;
    #pragma omp critical (input)
    {
      getline(in,ref);
    }
    if(exists_test(ref)){
      query_file_whole(ref);
    }
  }
}


void pVector(std::vector <gid> const &a) {
   std::cout << "The vector elements are : ";
   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';
}



void Index::toFile(const string &filename){
  //long totalSize;
  //DEBUG_MSG("Creating'" << filename << "'");
  ofstream outfile(filename);
  stringstream buffer;
  DEBUG_MSG("Index informations : " << endl
                << "- " << genome_numbers << " Indexed Genomes" << endl
                << "- " << Index::Buckets << " Buckets SIZE");

  pVector(*Buckets);
  //buffer << "### Index informations" << endl;
  //char* buffer = new char[totalSize];
  // write to outfile
  //DEBUG_MSG("Adding '" << buffer.str() << "' (" << buffer.str().length() << ")");
  buffer.str("");
  outfile.write(buffer.str().c_str(), buffer.str().length());
  // release dynamically-allocated memory
  //delete[] buffer;
  outfile.close();
}



atomic<uint32_t> genomes_downloaded(0);
atomic<uint64_t> bases_downloaded(0);


string exec(const char* cmd) {
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
	//~ cout<<str<<endl;
	uint lastposition(0);
	for(uint i(0);i<str.size()-3;++i){
		if(str[i]=='/'){
			lastposition=i;
		}
	}
	lastposition++;
	//~ cout<<str.substr(lastposition,str.size()-lastposition)<<endl;
	return str.substr(lastposition,str.size()-lastposition);
}


bool Index::Download_NCBI(const string& str, vector<uint64_t>& hashes){
	string cmd("wget -qO - "+str+"/"+get_name_ncbi(str)+"_genomic.fna.gz | gzip -d -c -");
	//~ cout<<cmd<<endl;
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
				//~ cout<<sequence.substr(0,10)<<endl;
				//~ cout<<sequence.substr(sequence.size()-10)<<endl;
				compute_sketch_kmer(sequence,hashes);
				bases_downloaded+=sequence.size();
				sequence.clear();
			}
		}
		result.clear();
    }
    if(something_to_eat){
		uint gtemp=++genomes_downloaded;
		if(gtemp%1000==0){
			cout<<"genomes ddl: "<<intToString(genomes_downloaded)<<" bases ddl: "<<intToString(bases_downloaded)<<endl;
		}
	}
	return something_to_eat;
}



//Hacky function to download genome from ncbi
void Index::Download_NCBI_fof(const string& fofncbi,const string& outfile){
	cout<<"go ncbi"<<endl;
	zstr::ifstream in(fofncbi);
	ofstream out(outfile);
	#pragma omp parallel
	while(not in.eof()){
		string ref;
		#pragma omp critical (in)
		{
			getline(in,ref);
		}
		if(ref.size()>5){
			vector<uint64_t> hashes(F,-1);
			if(Download_NCBI(ref,hashes)){
				//~ cout<<"yes"<<endl;
				#pragma omp critical (out)
				{
					out.write(reinterpret_cast<const char*>(&hashes[0]),F*8);
				}
			}
		}
		ref.clear();
	}
	cout<<"genomes ddl: "<<intToString(genomes_downloaded)<<" bases ddl: "<<intToString(bases_downloaded)<<endl;
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


