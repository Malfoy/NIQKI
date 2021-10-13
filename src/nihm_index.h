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
#include <unordered_set>
#include <sys/types.h>
#include <sys/stat.h>



using namespace std;
using gid = uint32_t;
using kmer = uint64_t;
using query_output=vector<pair<gid, uint32_t>>;



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



		Index();

		uint64_t nuc2int(char c) const;

		string kmer2str(uint64_t num,uint k) const ;

		uint64_t asm_log2(const uint64_t x) const;

		uint64_t nuc2intrc(char c) const;

		void update_kmer(kmer& min, char nuc)const;

		void update_kmer_RC(kmer& min, char nuc)const;

		kmer rcb(kmer min)const;

		void print_bin(uint64_t n,uint bits_to_print=64)const;

		kmer str2numstrand(const string& str)const;

		uint32_t get_fingerprint(uint64_t hashed)const;

		uint64_t revhash64 ( uint64_t x ) const;

		void compute_sketch(const string& reference, vector<int32_t>& sketch) const;

		//HERE we only select the minimal hashes without computing the HMH fingerprint
		void compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch) const;

		void insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id);

		//TODO Optmiser : compter les genomes id
		query_output query_sketch(const vector<int32_t>& sketch,uint32_t min_score=1) const;

		query_output query_sequence(const string& str,uint32_t min_score=1)const;

		query_output query_sequence(const string& str,double fraction)const;

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

		//HERE all the files of the fof are inserted as a separate entry in the index
		void insert_file_of_file_whole(const string& filestr);

		//HERE all the kmer of the file are put in a single sketch and Queried
		void query_file_whole(const string& filestr);

		void query_file_of_file_whole(const string& filestr);


};



#endif
