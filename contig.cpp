#ifndef CONGTIG_H
#define CONGTIG_H

#include <vector>
#include "read.cpp"
#include "settings.cpp"
#include "utils.cpp"

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
        int length = strlen(rhs._seq) +1; //+1 for \0 char
        _id = rhs._id;
        _seq = (char*)malloc(length);
        check_ptr(_seq);

        strncpy(_seq, rhs._seq, length);
        _qual = (char*)malloc(length);
        check_ptr(_qual);

        strncpy(_qual, rhs._qual, length);
    }

    //add new bases to end of sequence.  Give them quality of 1
    void append(char* added_seq){
        int added_length = strlen(added_seq);
        int current_length = size(); //need to store this before adding seq for changing quals
        int total_length = current_length + added_length;

        _seq = (char*)realloc(_seq, total_length+1);
        check_ptr(_seq);
        memcpy(&_seq[current_length], added_seq, added_length+1);

        _qual = (char*)realloc(_qual, total_length+1);
        check_ptr(_qual);
        for(int i=current_length; i<total_length; ++i){
            _qual[i] = '!';
        }
        _qual[total_length] = '\0';
    }

    unsigned int id() const{
        return _id;
    }

    void inc_qual(int pos){
        if(_qual[pos] < MAX_QUAL) {
            _qual[pos] = _qual[pos] + 1;
        }
    }

    //add more bases from front of sequence.  Give new bases quality of 1
    void prepend(char* added_seq){
        int added_length = strlen(added_seq);
        int current_length = size(); //need to grab this before adding quals
        int total_length = size() + added_length;

        _seq = (char*)realloc(_seq, total_length+1);
        check_ptr(_seq);

        memmove(&_seq[added_length], &_seq[0], current_length+1);
        memcpy(&_seq[0], added_seq, added_length);

        _qual = (char*)realloc(_qual, total_length+1);
        check_ptr(_qual);

        memmove(&_qual[added_length], &_qual[0], current_length+1);
        for(int i=0; i<added_length; ++i){
            _qual[i] = '!';
        }
    }

    char* qual() const {
        return _qual;
    }

    //compute and return the reverse complilment of the sequence
    //caller should ensure that returned pointer is freed.
    Contig* rev_comp(){
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
        check_ptr(rev_qual);

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
    
    char* seq() const{
        return _seq;
    }

    void set_id(unsigned int new_id){
        _id = new_id;
    }

    //should only be called by copy constructors & assignment operators.
    //usually quality should be itialized to 1 when new sequence is set.
    void set_qual(const char* new_qual){
        int length = strlen(new_qual);
        if(_qual == NULL){
            _qual = (char*)malloc(length+1);
        } else {
            _qual = (char*)realloc(_qual, length+1);
        }
        check_ptr(_qual);

        memcpy(_qual, new_qual, length+1);
    }

    //assign a new sequence.  Except in the case of copying from a different contig, 
    //it will usually be preferable to initialize the quality to 1 for all bases
    void set_seq(const char* new_seq, bool initialize_qual = true){
        int length = strlen(new_seq);

        if(_seq == NULL){
            _seq = (char*)malloc(length+1);
        } else {
            _seq = (char*)realloc(_seq, length+1);
        }
        check_ptr(_seq);

        memcpy(_seq, new_seq, length+1);

        if( initialize_qual ){
            if(_qual == NULL){
                _qual = (char*)malloc(length+1);
            } else {
                _qual = (char*)realloc(_qual, length+1);
            }
            check_ptr(_qual);

            for(int i=0; i<length; ++i){
                _qual[i] = 1 + QUAL_OFFSET;
            }
            _qual[length] = '\0';
        }
    }

    //replaces sequence and quality by those within the region sepecified.
    void set_substr(int start, int length){
        char* new_seq = (char*)malloc(length + 1);
        check_ptr(new_seq);
        char* new_qual = (char*)malloc(length + 1);
        check_ptr(new_qual);

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

    //return part fo the sequence.  caller must be sure to free pointer returned
    char* substr(int start, int length){
        char* sub = (char*)malloc(length+1);
        check_ptr(sub);
        memcpy(sub, &_seq[start], length);
        sub[length] = '\0';
        return sub;
    }

    //trim low quality bases from both sides of contig. shift & trim all aligned
    //reads if applicable.  If the contig has no bases above quality, unassemble
    //all of it's reads, as it will be deleted by the caller.
    void trim(char min_quality, vector<Read> &reads){
        min_quality += QUAL_OFFSET;

        //compute amount to trim on left side
        unsigned int start_pos = 0;
        while(start_pos < size() && _qual[start_pos] < min_quality){
            ++start_pos;
        }

        //compute amount to trim on right side
        unsigned int end_pos = size();
        while(end_pos > 0 && _qual[end_pos-1] < min_quality){
            --end_pos;
        }

        int new_length = end_pos - start_pos;

        //this contig is low quality and will be erased.  unassemble all of it's reads.
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

    //overloaded assignment operator for a deep copy.
    Contig& operator = (Contig &rhs){
        _id = rhs._id;
        int length = rhs.size() +1;

        _seq = (char*)malloc(length);
        check_ptr(_seq);
        strncpy(_seq, rhs._seq, length);

        _qual = (char*)malloc(length);
        check_ptr(_qual);
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
