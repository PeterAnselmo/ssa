#ifndef CONGTIG_H
#define CONGTIG_H

#include <algorithm>
#include <vector>
#include "read.cpp"

using namespace std;

class Contig {
private:
    char* _seq;
    char* _qual;
    unsigned int _id;

public:

    Contig(){
        _seq = NULL;
        _qual = NULL;
    }
    Contig(int id){
        _id = id;
        _seq = NULL;
        _qual = NULL;
    }
    Contig(int id, char* new_seq){
        _id = id;
        _seq = NULL;
        _qual = NULL;
        set_seq(new_seq);
    }
    Contig(const Contig &rhs){
        int length = strlen(rhs._seq) +1;
        _id = rhs._id;
        _seq = (char*)malloc(length);
        strncpy(_seq, rhs._seq, length);
        _qual = (char*)malloc(length);
        strncpy(_qual, rhs._qual, length);
    }

    void append(char* added_seq){
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
    }

    unsigned int id() const{
        return _id;
    }

    void inc_qual(int pos){
        if(_qual[pos] < 126) {
            _qual[pos] = _qual[pos] + 1;
        }
    }

    void prepend(char* added_seq){
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
    }

    char* qual() const {
        return _qual;
    }
    
    char* seq() const{
        return _seq;
    }

    void set_seq(const char* new_seq){
        int length = strlen(new_seq);
        if(_seq == NULL){
            _seq = (char*)malloc(length+1);
        } else {
            _seq = (char*)realloc(_seq, length+1);
        }
        if(_seq == NULL){
            printf("Memory Allocation Error.\n");
            exit(1);
        }
        memcpy(_seq, new_seq, length+1);
        if(_qual == NULL){
            _qual = (char*)malloc(length+1);
        } else {
            _qual = (char*)realloc(_qual, length+1);
        }
        for(int i=0; i<length; ++i){
            _qual[i] = '!';
        }
        _qual[length] = '\0';
    }

    void set_substr(int start, int length){
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
    void shift_aligned_reads(unsigned int distance, vector<Read> &reads){
        vector<Read>::iterator read;
        for(read = reads.begin(); read != reads.end(); ++read){
            if(read->assembled() && read->contig() == _id){
                int overlap = distance - read->position();
                if( overlap > 0){
                    char *gapped_substr = read->gapped_substr(overlap);
                    read->set_gapped_seq(gapped_substr);
                    free(gapped_substr);
                    read->set_position(0);
                } else {
                    read->set_position(read->position() - distance);
                }
            }
        }
    }

    unsigned int size(){
        return strlen(_seq);
    }

    char* substr(int start, int length){
        char* sub = (char*)malloc(length+1);
        memcpy(sub, &_seq[start], length);
        sub[length] = '\0';
        return sub;
    }

    void trim(char min_quality, vector<Read> &reads){
        min_quality += 33;

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
            vector<Read>::iterator read;
            for(read = reads.begin(); read != reads.end(); ++read){
                if( read->assembled() && read->contig() == _id ){
                    read->unassemble();
                }
            }
        } else {
            if(start_pos > 0){
                shift_aligned_reads(start_pos, reads);
            }
        }
        set_substr(start_pos, new_length);
    }

    //if more was added to the start of the contig, it will be necessarry
    //to move all aligned reads to the right
    void unshift_aligned_reads(unsigned int distance, vector<Read> &reads){
        vector<Read>::iterator read;
        for(read = reads.begin(); read != reads.end(); ++read){
            if(read->assembled() && read->contig() == _id){
                read->set_position(read->position() + distance);
            }
        }
    }

    Contig& operator = (Contig &rhs){
        _id = rhs._id;
        int length = rhs.size() +1;
        _seq = (char*)malloc(length);
        strncpy(_seq, rhs._seq, length);
        _qual = (char*)malloc(length);
        strncpy(_seq, rhs._seq, length);
        return *this;
    }

    ~Contig(){
        if(_seq != NULL){
            free(_seq);
        }
        if(_qual != NULL){
            free(_qual);
        }
    }

};

#endif
