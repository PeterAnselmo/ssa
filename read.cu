#ifndef READ_CU
#define READ_CU

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include "cudastring.cu"

class Read {
private:
    
    char* _description;
    char* _seq;
    char* _qual;
    char* _gapped_seq;

    bool _assembled;
    unsigned int _contig;
    unsigned int _position;
public:

    __host__ __device__ Read(){
        _description = NULL;
        _seq = NULL;
        _qual = NULL;
        _gapped_seq = NULL;

        _assembled = false;
        _contig = 0;
        _position = 0;
    }
    __host__ __device__ Read(const Read &rhs){
        _description = NULL;
        _seq = NULL;
        _qual = NULL;
        _gapped_seq = NULL;

        set_description(rhs._description);
        set_seq(rhs._seq);
        set_qual(rhs._qual);
        set_gapped_seq(rhs._gapped_seq);

        _assembled = rhs._assembled;
        _contig = rhs._contig;
        _position = rhs._position;
    }

    __host__ __device__ bool assembled() const{
        return _assembled;
    }

    __host__ __device__ void assemble(unsigned int contig_id, unsigned int pos){
        _assembled = true;
        _contig = contig_id;
        _position = pos;
        set_gapped_seq(_seq);
    }

    __host__ __device__ unsigned int contig() const {
        return _contig;
    }

    __host__ __device__ const char* description() const{
        return _description;
    }

    __host__ __device__ unsigned int find(char search){
        #ifndef __CUDA_ARCH__
            for(unsigned int i=0; i<strlen(_seq); ++i){
                if(_seq[i] == search){
                    return i;
                }
            }
            return size();
        #else
            for(unsigned int i=0; i<cudaStrlen(_seq); ++i){
                if(_seq[i] == search){
                    return i;
                }
            }
            return size();
        #endif
    }

    __host__ __device__ char* gapped_seq(){
        return _gapped_seq;
    }

    __host__ __device__ unsigned int gapped_size() const {
        return strlen(_gapped_seq);
    }

    __host__ __device__ char* gapped_substr(int pos, int length) const {
        char *sub = (char *)malloc(sizeof(char)*(length+1));
        memcpy(sub, &_gapped_seq[pos], length);
        sub[length] = '\0';
        return sub;
    }
    __host__ __device__ char* gapped_substr(int pos) const {
        int length = size() - pos;
        return gapped_substr(pos, length);
    }

    __host__ __device__ unsigned int position(){
        return _position;
    }

    __host__ __device__ const char* qual() const{
        return _qual;
    }

    __host__ __device__ char* rev_substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        char *rev = rev_comp();
        memcpy(sub, &rev[pos], length);
        free(rev);
        sub[length] = '\0';
        return sub;
    }

    __host__ __device__ char* rev_comp() const{
        char *rev = (char *)malloc(sizeof(char)*(size()+1));

        for(unsigned int i=0; i <size(); ++i){
            int pos = size()-1-i;
            switch(_seq[pos]){
                case 'A':
                    rev[i] = 'T';
                break;
                case 'T':
                    rev[i] = 'A';
                break;
                case 'C':
                    rev[i] = 'G';
                break;
                case 'G':
                    rev[i] = 'C';
                break;
                case 'N':
                    rev[i] = 'N';
                break;
            }
        }
        rev[size()] = '\0';
        return rev;
    }

    __host__ __device__ const char* seq() const{
        return _seq;
    }

    __host__ __device__ void set_rev_comp(){
        char* rev = rev_comp();
        set_seq(rev);
        free(rev);
        //todo - dry up string reversal
        //reverse(_qual.begin(), _qual.end());
    }

    __host__ __device__ void set_description(const char *new_description){
        #ifndef __CUDA_ARCH__
            if(_description == NULL){
                _description = (char*)malloc(strlen(new_description)+1);
            } else {
                _description = (char*)realloc(_description, strlen(new_description)+1);
            }
            if(_description == NULL){
                fprintf(stderr, "Memory Allocation Error\n");
                exit(1);
            }
            strncpy(_description, new_description, strlen(new_description)+1);
        #else
            if(_description == NULL){
                _description = (char*)malloc(cudaStrlen(new_description)+1);
            } else {
                free(_description);
                _description = (char*)malloc(cudaStrlen(new_description)+1);
            }
            cudaStrncpy(_description, new_description, cudaStrlen(new_description)+1);
        #endif
    }

    __host__ __device__ void set_gapped_seq(const char* new_gapped_seq){
        #ifndef __CUDA_ARCH__
            if(_gapped_seq == NULL){
                _gapped_seq = (char*)malloc(strlen(new_gapped_seq)+1);
            } else {
                _gapped_seq = (char*)realloc(_gapped_seq, strlen(new_gapped_seq)+1);
            }
            strncpy(_gapped_seq, new_gapped_seq, strlen(new_gapped_seq)+1);
        #else
            if(_gapped_seq == NULL){
                _gapped_seq = (char*)malloc(cudaStrlen(new_gapped_seq)+1);
            } else {
                free(_gapped_seq);
                _gapped_seq = (char*)malloc(cudaStrlen(new_gapped_seq)+1);
            }
            cudaStrncpy(_gapped_seq, new_gapped_seq, cudaStrlen(new_gapped_seq)+1);
        #endif
    }

    __host__ __device__ void set_position(unsigned int new_pos){
        _position = new_pos;
    }

    __host__ __device__ void set_qual(const char* new_qual){
        #ifndef __CUDA_ARCH__
            if(_qual == NULL){
                _qual = (char*)malloc(strlen(new_qual)+1);
            } else {
                _qual = (char*)realloc(_qual, strlen(new_qual)+1);
            }
            strncpy(_qual, new_qual, strlen(new_qual)+1);
        #else
            if(_qual == NULL){
                _qual = (char*)malloc(cudaStrlen(new_qual)+1);
            } else {
                free(_qual);
                _qual = (char*)malloc(cudaStrlen(new_qual)+1);
            }
            cudaStrncpy(_qual, new_qual, cudaStrlen(new_qual)+1);
        #endif
    }

    __host__ __device__ void set_seq(const char* new_seq){
        #ifndef __CUDA_ARCH__
            if(_seq == NULL){
                _seq = (char*)malloc(strlen(new_seq)+1);
            } else {
                _seq = (char*)realloc(_seq, strlen(new_seq)+1);
            }
            strncpy(_seq, new_seq, strlen(new_seq)+1);
        #else
            if(_seq == NULL){
                _seq = (char*)malloc(cudaStrlen(new_seq)+1);
            } else {
                free(_seq);
                _seq = (char*)malloc(cudaStrlen(new_seq)+1);
            }
            cudaStrncpy(_seq, new_seq, cudaStrlen(new_seq)+1);
        #endif
    }

    __host__ __device__ unsigned int size() const {
        #ifndef __CUDA_ARCH__
            return strlen(_seq);
        #else
            return cudaStrlen(_seq);
        #endif
    }

    __host__ __device__ char* substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        memcpy(sub, &_seq[pos], length);
        sub[length] = '\0';
        return sub;
    }

    __host__ __device__ char* substr(int pos) const {
        int length = size()-pos;
        return substr(pos, length);
    }

    __host__ __device__ void trim(int num_bases = 2){
        char *new_seq;

        #ifndef __CUDA_ARCH__
            int new_size = strlen(_seq) - 2*num_bases;
        #else
            int new_size = cudaStrlen(_seq) - 2*num_bases;
        #endif

        //read is too short to trim desired amount
        if( new_size <= 0){
            return;
        }

        new_seq = (char*)malloc(new_size+1);
        for(int i=0; i<new_size; ++i){
            new_seq[i] = _seq[i+num_bases];
        }
        new_seq[new_size] = '\0';
        free(_seq);
        _seq = new_seq;
    }

    __host__ __device__ void unassemble(){
        _assembled = false;
    }

    bool operator < (const Read &rhs) const {

        if(_assembled){
            if(rhs._assembled){
                if(_contig == rhs._contig){
                    return _position < rhs._position;
                }
                return _contig < rhs._contig;
            }
            return true;
        }
        return false;
    }

    Read& operator = (const Read &rhs){
        free_strings();

        set_description(rhs._description);
        set_seq(rhs._seq);
        set_qual(rhs._qual);
        set_gapped_seq(rhs._gapped_seq);

        _assembled = rhs._assembled;
        _contig = rhs._contig;
        _position = rhs._position;

        return *this;
    }

    __host__ __device__ ~Read(){
        free_strings();
    }

private:
    __host__ __device__ void free_strings(){
        if(_description != NULL){
            free(_description);
            _description = NULL;
        }
        if(_seq != NULL){
            free(_seq);
            _seq = NULL;
        }
        if(_qual != NULL){
            free(_qual);
            _qual = NULL;
        }
        if(_gapped_seq != NULL){
            free(_gapped_seq);
            _gapped_seq = NULL;
        }
    }
};

#endif
