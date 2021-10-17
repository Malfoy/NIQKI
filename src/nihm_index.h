#ifndef __NIHM_INDEX_H__
#define __NIHM_INDEX_H__



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <algorithm>
#include <mutex>
#include <unordered_set>
#include <sys/types.h>
#include <sys/stat.h>
#include <omp.h>



using namespace std;
using gid = uint32_t;
using kmer = uint64_t;
using query_output=vector<pair<gid, uint32_t>>;
const uint32_t mutex_number=65536;


class Index {
  public:
    //CONSTANTS
    uint32_t K;//kmer size
    uint32_t F;//fingerprint used
    uint32_t W;//fingerprint size
    uint32_t M;//Minhash size
    uint32_t mask_M;//minhash mask
    uint32_t H;//Hll size (M+H=W)
    uint32_t maximal_remainder;//2^H-1
    uint32_t lF;//log2(F)
    uint32_t fingerprint_range;//2^w
    uint64_t mask_fingerprint;//2^(64-lf)-1
    uint64_t expected_gemome_size;
    uint64_t offsetUpdatekmer;
    //VARIABLE
    uint32_t genome_numbers;//Number of genomes
    vector<gid>* Buckets;
    omp_lock_t lock[mutex_number];



    /**
     * \brief Default constructor.
     */
    Index();
    Index(uint32_t lF, uint32_t K, uint32_t W, uint32_t H);

    /**
     * \brief Destructor.
     */
    ~Index();

    uint64_t nuc2int(char c) const;

    string kmer2str(uint64_t num,uint k) const ;

    uint64_t asm_log2(const uint64_t x) const;

    uint64_t nuc2intrc(char c) const;

    void update_kmer(kmer& min, char nuc)const;

    void update_kmer_RC(kmer& min, char nuc)const;

    kmer rcb(kmer min)const;

    kmer str2numstrand(const string& str)const;

    uint32_t get_fingerprint(uint64_t hashed)const;

    uint64_t revhash64 ( uint64_t x ) const;

    void compute_sketch(const string& reference, vector<int32_t>& sketch) const;

    //HERE we only select the minimal hashes without computing the HMH fingerprint
    void compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch) const;

    void insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id);

    inline void print_bin(uint64_t n,uint bits_to_print=64) const{
        uint64_t mask=1;
        mask<<=bits_to_print-1;
        for(uint i(0);i<bits_to_print;++i){
            cout<<n/mask;
            if(n/mask==1){n-=mask;}
            mask>>=1;
        }
        cout<<"\n";
    }

    inline query_output query_sketch(const vector<int32_t>& sketch,uint32_t min_score=1)const {
        query_output result;
        unordered_map<gid,uint32_t> counts;
        for(uint i(0);i<F;++i){
            if(sketch[i]<fingerprint_range and sketch[i]>0){
                for(uint j(0);j<Buckets[sketch[i]+i*fingerprint_range].size();++j){
                    counts[Buckets[sketch[i]+i*fingerprint_range][j]]++;
                }
            }
        }
        for(auto it=counts.begin();it!=counts.end();it++) {
            if(it->second>=min_score){
                result.push_back({it->second,it->first});
            }
        }
        sort(result.begin(),result.end(),greater<pair<uint32_t,uint32_t>>());
        return result;
    }



    inline query_output query_sequence(const string& str,uint32_t min_score=1)const {
        vector<int32_t> sketch;
        compute_sketch(str,sketch);
        return query_sketch(sketch,min_score);
    }



    inline query_output query_sequence(const string& str,double fraction)const {
        return query_sequence(str,F*fraction);
    }


    void insert_sequence(const string& str,uint32_t genome_id);

    //all the lines from the file is considered as a separate entry with a different identifier
    //TODO HANDLE FASTQ multiFASTA
    void insert_file_lines(const string& filestr);

    void query_file_lines(const string& filestr);

    void merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2);

    inline bool exists_test (const std::string& name) {
      struct stat buffer;
      return (stat (name.c_str(), &buffer) == 0);
    }



    //HERE all the kmer of the file are put in a single sketch and inserted
    void insert_file_whole(const string& filestr);

    void insert_file_whole(const string& filestr,uint32_t identifier);

    //HERE all the files of the fof are inserted as a separate entry in the index
    void insert_file_of_file_whole(const string& filestr);

    //HERE all the kmer of the file are put in a single sketch and Queried
    void query_file_whole(const string& filestr);

    void query_file_of_file_whole(const string& filestr);

     /**
     * \brief Write the Index object into the given file.
     *
     * \param filename The file to write.
     */
    void toFile(const string &filename);

};



#endif
