#ifndef CONGTIG_H
#define CONGTIG_H

#include "read.cu"
#include "settings.cu"

using namespace std;

class Contig {
private:
    char* _seq;
    char* _qual;
    unsigned int _id;

public:

    __host__ __device__ Contig(){
        _seq = NULL;
        _qual = NULL;
    }
    __host__ __device__ Contig(int id){
        _id = id;
        _seq = NULL;
        _qual = NULL;
    }
    __host__ __device__ Contig(int id, char* new_seq){
        _id = id;
        _seq = NULL;
        _qual = NULL;
        set_seq(new_seq);
    }
    __host__ __device__ Contig(const Contig &rhs){
        int length = strlen(rhs._seq) +1;
        _id = rhs._id;
        _seq = (char*)malloc(length);
        strncpy(_seq, rhs._seq, length);
        _qual = (char*)malloc(length);
        strncpy(_qual, rhs._qual, length);
    }

    __host__ __device__ void append(char* added_seq){
        #ifndef __CUDA_ARCH__
            int added_length = strlen(added_seq);
            int current_length = size(); //need to grab this before changing quals
            int total_length = current_length + added_length;

            _seq = (char*)realloc(_seq, total_length+1);
            memcpy(&_seq[current_length], added_seq, added_length+1);
            //_seq[total_length] = '\0';

            _qual = (char*)realloc(_qual, total_length+1);
            for(int i=current_length; i<total_length; ++i){
                _qual[i] = '!';
            }
            _qual[total_length] = '\0';
        #else
            int added_length = cudaStrlen(added_seq);
            int current_length = size(); //need to grab this before changing quals
            int total_length = current_length + added_length;

            free(_seq);
            _seq = (char*)malloc(total_length+1);
            memcpy(&_seq[current_length], added_seq, added_length+1);
            //_seq[total_length] = '\0';

            free(_qual);
            _qual = (char*)malloc(total_length+1);
            for(int i=current_length; i<total_length; ++i){
                _qual[i] = '!';
            }
            _qual[total_length] = '\0';

        #endif
    }

    __host__ __device__ unsigned int id() const{
        return _id;
    }

    __host__ __device__ void inc_qual(int pos){
        if(_qual[pos] < 126) {
            _qual[pos] = _qual[pos] + 1;
        }
    }

    __host__ __device__ void prepend(char* added_seq){
        #ifndef __CUDA_ARCH__
            int added_length = strlen(added_seq);
            int current_length = size(); //need to grab this before changing quals
            int total_length = size() + added_length;

            _seq = (char*)realloc(_seq, total_length+1);
            if(_seq == NULL){
                printf("Memory Allocation Error.\n");
                exit(1);
            }
            memmove(&_seq[added_length], &_seq[0], current_length+1);
            memcpy(&_seq[0], added_seq, added_length);

            _qual = (char*)realloc(_qual, total_length+1);
            memmove(&_qual[added_length], &_qual[0], current_length+1);
            for(int i=0; i<added_length; ++i){
                _qual[i] = '!';
            }
        #else
            int added_length = cudaStrlen(added_seq);
            int current_length = size(); //need to grab this before changing quals
            int total_length = size() + added_length;

            char* new_seq = (char*)malloc(total_length+1);
            for(int i=0; i<added_length; ++i){
                new_seq[i] = added_seq[i];
            }
            for(int i=0; i<=current_length; ++i){
                new_seq[i+added_length] = _seq[i];
            }

            char* new_qual = (char*)malloc(total_length+1);
            for(int i=0; i<added_length; ++i){
                _qual[i] = '!';
            }
            for(int i=0; i<=current_length; ++i){
                new_qual[i+added_length] = _qual[i];
            }
            free(_seq);
            free(_qual);
            _seq = new_seq;
            _qual = new_qual;
        #endif
    }

    __host__ __device__ char* qual() const {
        return _qual;
    }

    __host__ __device__ Contig* rev_comp(){
        char *rev = (char *)malloc(size()+1);

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

        char *rev_qual = (char*)malloc(size()+1);
        for(unsigned int i=0; i<size(); ++i){
            int pos = size()-1-i;
            rev_qual[i] = _qual[pos];
        }
        rev_qual[size()] = '\0';

        //Contig* temp = (Contig*)malloc(sizeof(Contig));
        Contig* temp = new Contig();
        temp->set_id(_id);
        temp->set_seq(rev, false);
        temp->set_qual(rev_qual);

        free(rev);
        free(rev_qual);

        return temp;
    }
    
    __host__ __device__ char* seq() const{
        return _seq;
    }

    __host__ __device__ void set_id(unsigned int new_id){
        _id = new_id;
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

    __host__ __device__ void set_seq(const char* new_seq, bool initialize_qual = true){
        #ifndef __CUDA_ARCH__
            int length = strlen(new_seq);
            if(_seq == NULL){
                _seq = (char*)malloc(length+1);
            } else {
                _seq = (char*)realloc(_seq, length+1);
            }
            strncpy(_seq, new_seq, length+1);
        #else
            int length = cudaStrlen(new_seq);
            if(_seq == NULL){
                _seq = (char*)malloc(length+1);
            } else {
                free(_seq);
                _seq = (char*)malloc(length+1);
            }
            cudaStrncpy(_seq, new_seq, length+1);
        #endif
        if( initialize_qual ){
            #ifndef __CUDA_ARCH__
                if(_qual == NULL){
                    _qual = (char*)malloc(length+1);
                } else {
                    _qual = (char*)realloc(_qual, length+1);
                }
            #else
                if(_qual == NULL){
                    _qual = (char*)malloc(length+1);
                } else {
                    free(_qual);
                    _qual = (char*)malloc(length+1);
                }
            #endif
            for(int i=0; i<length; ++i){
                _qual[i] = '!';
            }
            _qual[length] = '\0';
        }
    }

    __host__ __device__ void set_substr(int start, int length){
        char* new_seq = (char*)malloc(length + 1);
        char* new_qual = (char*)malloc(length + 1);

        memcpy(new_seq, &_seq[start], length);
        memcpy(new_qual, &_qual[start], length);
        new_seq[length] = '\0';
        new_qual[length] = '\0';

        free(_seq);
        free(_qual);

        _seq = new_seq;
        _qual = new_qual;

    }
    
    //if sequence was removed from the start of the contig, it will be necessary
    //to move all reads to the left, and possibly trim them
    __host__ __device__ void shift_aligned_reads(unsigned int distance, Read *reads, int num_reads){
        for(int i=0; i<num_reads; ++i){
            Read read = reads[i];
            if(read.assembled() && read.contig() == _id){
                int overlap = distance - read.position();
                if( overlap > 0){
                    char *gapped_substr = read.gapped_substr(overlap);
                    read.set_gapped_seq(gapped_substr);
                    free(gapped_substr);
                    read.set_position(0);
                } else {
                    read.set_position(read.position() - distance);
                }
            }
        }
    }

    __host__ __device__ unsigned int size(){
        #ifndef __CUDA_ARCH__
            return strlen(_seq);
        #else
            return cudaStrlen(_seq);
        #endif
    }

    __host__ __device__ char* substr(int start, int length){
        char* sub = (char*)malloc(length+1);
        memcpy(sub, &_seq[start], length);
        sub[length] = '\0';
        return sub;
    }

    __host__ __device__ void trim(char min_quality, Read reads[], int num_reads){
        min_quality += QUAL_OFFSET;

        //trim_left_size
        unsigned int start_pos = 0;
        while(start_pos < size() && _qual[start_pos] < min_quality){
            ++start_pos;
        }

        //trim_right_size
        unsigned int end_pos = size();
        while(end_pos > 0 && _qual[end_pos-1] < min_quality){
            --end_pos;
        }

        int new_length = end_pos - start_pos;
        if(new_length <= 0){ 
            new_length = 0; 
            for(int i=0; i<num_reads; ++i){
                if( reads[i].assembled() && reads[i].contig() == _id ){
                    reads[i].unassemble();
                }
            }
        } else {
            if(start_pos > 0){
                shift_aligned_reads(start_pos, reads, num_reads);
            }
        }
        set_substr(start_pos, new_length);
    }

    //if more was added to the start of the contig, it will be necessarry
    //to move all aligned reads to the right
    __host__ __device__ void unshift_aligned_reads(unsigned int distance, Read* reads, int num_reads){
        for(int i=0; i<num_reads; ++i){
            if(reads[i].assembled() && reads[i].contig() == _id){
                reads[i].set_position(reads[i].position() + distance);
            }
        }
    }

    __host__ __device__ Contig& operator = (Contig &rhs){
        _id = rhs._id;
        int length = rhs.size() +1;
        _seq = (char*)malloc(length);
        strncpy(_seq, rhs._seq, length);
        _qual = (char*)malloc(length);
        strncpy(_seq, rhs._seq, length);
        return *this;
    }

    __host__ __device__ ~Contig(){
        if(_seq != NULL){
            free(_seq);
        }
        if(_qual != NULL){
            free(_qual);
        }
    }

};

#endif
