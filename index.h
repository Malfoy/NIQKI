#ifndef index
#define index



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



class Index{
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



    Index(uint32_t ilF=10, uint32_t iK=31,uint32_t iW=8,uint32_t iH=4){
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
    }



    ~Index(){
        delete[] Buckets;
    }




    uint64_t nuc2int(char c)const{
        switch(c){
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return 0;
        }
        exit(0);
        return 0;
    }



    string kmer2str(uint64_t num,uint k)const{
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




    uint64_t asm_log2(const uint64_t x) const{
        uint64_t y;
        asm ( "\tbsr %1, %0\n"
            : "=r"(y)
            : "r" (x)
        );
        return y;
    }


    uint64_t mylog2 (uint64_t val) const{
        if (val == 0) return 0;
        if (val == 1) return 0;
        uint64_t ret = 0;
        while (val > 1) {
            val >>= 1;
            ret++;
        }
        return ret;
    }



    uint64_t nuc2intrc(char c)const{
        switch(c){
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




    void update_kmer(kmer& min, char nuc)const{
        min<<=2;
        min+=nuc2int(nuc);
        min%=offsetUpdatekmer;
    }



    void update_kmer_RC(kmer& min, char nuc)const{
        min>>=2;
        min+=(nuc2intrc(nuc)<<(2*K-2));
    }



    kmer rcb(kmer min)const{
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


    void print_bin(uint64_t n,uint bits_to_print=64)const{
        uint64_t mask=1;
        mask<<=bits_to_print-1;
        for(uint i(0);i<bits_to_print;++i){
            cout<<n/mask;
            if(n/mask==1){n-=mask;}
            mask>>=1;
        }
        cout<<"\n";
    }



    kmer str2numstrand(const string& str)const{
        uint64_t res(0);
        for(uint i(0);i<str.size();i++){
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



    uint32_t get_fingerprint(uint64_t hashed)const{
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
        //   cout<<"get_fingerprintOUT"<<endl;
        // cout<<remaining_nonzero<<endl;
        //  print_bin(result,8);
        //  cin.get();
        return result;
    }



    uint64_t revhash64 ( uint64_t x ) const {
        x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
        x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
        x = ( ( x >> 32 ) ^ x );
        return x;
    }



    void compute_sketch(const string& reference, vector<int32_t>& sketch)const{
        if(sketch.size()!=F){
            sketch.resize(F,-1);
        }
        kmer S_kmer(str2numstrand(reference.substr(0,K-1)));
	    kmer RC_kmer(rcb(S_kmer));
        for(uint i(0);i+K<reference.size();++i){
            update_kmer(S_kmer,reference[i+K-1]);
            update_kmer_RC(RC_kmer,reference[i+K-1]);
            kmer canon(min(S_kmer,RC_kmer));
            uint64_t hashed=revhash64(canon);
            uint32_t bucket_id(hashed>>(64-lF));
            // print_bin(hashed);
            hashed&=mask_fingerprint;
            // print_bin(hashed);
            // print_bin(sketch[bucket_id]);
            hashed=get_fingerprint(hashed);
            if(sketch[bucket_id]>hashed){
                sketch[bucket_id]=hashed;
            }
        }
    }



    //HERE we only select the minimal hashes without computing the HMH fingerprint
     void compute_sketch_kmer(const string& reference, vector<uint64_t>& sketch)const{
        // cout<<"compute_sketch_kmer"<<endl;
        if(sketch.size()!=F){
            sketch.resize(F,-1);
        }
        kmer S_kmer(str2numstrand(reference.substr(0,K-1)));
	    kmer RC_kmer(rcb(S_kmer));
        for(uint i(0);i+K<reference.size();++i){
            update_kmer(S_kmer,reference[i+K-1]);
            update_kmer_RC(RC_kmer,reference[i+K-1]);
            kmer canon(min(S_kmer,RC_kmer));
            uint64_t hashed=revhash64(canon);
            uint32_t bucket_id(hashed>>(64-lF));
            hashed&=mask_fingerprint;
            // cout<<bucket;_id<<" "<<sketch.size()<<endl;
            if(sketch[bucket_id]>hashed){
                sketch[bucket_id]=hashed;
            }
        }
        // cout<<"compute_sketch_kmer end"<<endl;
    }



    void insert_sketch(const vector<int32_t>& sketch,uint32_t genome_id){
        for(uint i(0);i<F;++i){
            if(sketch[i]<fingerprint_range and sketch[i]>=0){
                Buckets[sketch[i]+i*fingerprint_range].push_back(genome_id);
            }
        }
    }



    query_output query_sketch(const vector<int32_t>& sketch,uint32_t min_score=1)const {
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



    query_output query_sequence(const string& str,uint32_t min_score=1)const {
        vector<int32_t> sketch;
        compute_sketch(str,sketch);
        return query_sketch(sketch,min_score);
    }



    query_output query_sequence(const string& str,double fraction)const {
        return query_sequence(str,F*fraction);
    }



    void insert_sequence(const string& str,uint32_t genome_id){
        vector<int32_t> sketch;
        compute_sketch(str,sketch);
        insert_sketch(sketch,genome_id);
    }



    //all the lines from the file is considered as a separate entry with a different identifier
    //TODO HANDLE FASTQ multiFASTA
    void insert_file_lines(const string& filestr){
        ifstream in(filestr);
        string ref,head;
        while(not in.eof()){
			{
				getline(in,head);
				getline(in,ref);
			}
            if(ref.size()>K){
                insert_sequence(ref,genome_numbers);
                genome_numbers++;
            }
            ref.clear();
           
        }
    }



    void query_file_lines(const string& filestr){
        ifstream in(filestr);
        string ref,head;
        while(not in.eof()){
            {
                getline(in,head);
                getline(in,ref);
            }
            if(ref.size()>K){
                auto out(query_sequence(ref));
                ref.clear();
                cout<<"out: "<<endl;
                for(uint32_t i=0;i<min((uint)5,(uint)out.size());i++){
                    cout<<out[i].second<<" "<<out[i].first<<endl;
                }
            }

        }
    }


        
    void merge_sketch( vector<int32_t>& sketch1,const vector<int32_t>& sketch2){
        for(uint i(0);i<sketch1.size();++i){
            sketch1[i]=min(sketch1[i],sketch2[i]);
        }
    }



inline bool exists_test (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}



//HERE all the kmer of the file are put in a single sketch and inserted
    void insert_file_whole(const string& filestr){
        ifstream in(filestr);
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
            ref.clear();
        }   
        // cout<<"get fingerprint"<<endl;
        for(uint i(0);i<F;++i){
            sketch[i]=get_fingerprint(kmer_sketch[i]);
        }
        // cout<<"get fingerprintOK"<<endl;

        insert_sketch(sketch,genome_numbers);;
        genome_numbers++;
    }



//HERE all the files of the fof are inserted as a separate entry in the index
    void insert_file_of_file_whole(const string& filestr){
        ifstream in(filestr);
        string ref;

        while(not in.eof()){
			getline(in,ref);
            if(exists_test(ref)){
                insert_file_whole(ref);
            }
        }
    }



    //HERE all the kmer of the file are put in a single sketch and Queried
    void query_file_whole(const string& filestr){
        ifstream in(filestr);
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
        for(uint i(0);i<F;++i){
            // cout<<i<<" "<<sketch.size()<<endl;
            sketch[i]=get_fingerprint(kmer_sketch[i]);
        }
        // cout<<"la fingernbtoint"<<endl;
        auto out(query_sketch(sketch,10000));
        cout<<out.size()<<endl;
    }



    void query_file_of_file_whole(const string& filestr){
        ifstream in(filestr);
        string ref;

        while(not in.eof()){
            getline(in,ref);
            if(exists_test(ref)){
                query_file_whole(ref);
            }
        }
    }



};




#endif