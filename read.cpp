#ifndef READ_CU
#define READ_CU

#include <stdio.h>
#include <cstdlib>
#include <cstring>

class Read {
private:
    
    char* _description;
    char* _seq;
    char* _qual;
    bool _assembled;
    unsigned int _contig;
    unsigned int _position;
public:
    char* gapped_seq;

    Read(){
        _assembled = false;
        _description = NULL;
        _seq = NULL;
        _qual = NULL;
    }
    Read(const Read &rhs){
        _description = NULL;
        _seq = NULL;
        _qual = NULL;
        set_description(rhs._description);
        set_seq(rhs._seq);
        set_qual(rhs._qual);
        _assembled = rhs._assembled;
    }

    bool assembled() const{
        return _assembled;
    }

    void assemble(unsigned int contig_id, unsigned int pos){
        _assembled = true;
        _contig = contig_id;
        _position = pos;
        gapped_seq = _seq;
    }

    unsigned int contig() const {
        return _contig;
    }

    const char* description() const{
        return _description;
    }

    unsigned int find(char search){
        for(unsigned int i=0; i<strlen(_seq); ++i){
            if(_seq[i] == search){
                return i;
            }
        }
        return size();
    }

    unsigned int gapped_size() const {
        return strlen(gapped_seq);
    }

    char* gapped_substr(int pos, int length) const {
        char *sub = (char *)malloc(sizeof(char)*(length+1));
        memcpy(sub, &gapped_seq[pos], length);
        sub[length] = '\0';
        return sub;
    }
    char* gapped_substr(int pos) const {
        int length = size() - pos;
        return gapped_substr(pos, length);
    }

    unsigned int position(){
        return _position;
    }

    const char* qual() const{
        return _qual;
    }

    char* rev_substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        char *rev = rev_comp();
        memcpy(sub, &rev[pos], length);
        free(rev);
        sub[length] = '\0';
        return sub;
    }

    char* rev_comp() const{
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

    const char* seq() const{
        return _seq;
    }

    void set_rev_comp(){
        strcpy(_seq, rev_comp());
        //todo - dry up string reversal
        //reverse(_qual.begin(), _qual.end());
    }

    void set_description(const char *new_description){
        if(_description == NULL){
            _description = (char *)malloc(strlen(new_description)+1);
        } else {
            _description = (char *)realloc(_description, strlen(new_description)+1);
        }
        if(_description == NULL){
            fprintf(stderr, "Memory Allocation Error\n");
            exit(1);
        }
        strncpy(_description, new_description, strlen(new_description)+1);
    }

    void set_position(unsigned int new_pos){
        _position = new_pos;
    }

    void set_qual(const char* new_qual){
        if(_qual == NULL){
            _qual = (char *)malloc(strlen(new_qual)+1);
        } else {
            _qual = (char *)realloc(_qual, strlen(new_qual)+1);
        }

        if(_qual == NULL){
            fprintf(stderr, "Memory Allocation Error\n");
            exit(1);
        }
        strncpy(_qual, new_qual,strlen(new_qual)+1);
    }

    void set_seq(const char* new_seq){
        if(_seq == NULL){
            _seq = (char *)malloc(strlen(new_seq)+1);
        } else {
            _seq = (char *)realloc(_seq, strlen(new_seq)+1);
        }
        if(_seq == NULL){
            fprintf(stderr, "Memory Allocation Error\n");
            exit(1);
        }
        strncpy(_seq, new_seq, strlen(new_seq)+1);
    }

    unsigned int size() const {
        return strlen(_seq);
    }

    char* substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        memcpy(sub, &_seq[pos], length);
        sub[length] = '\0';
        return sub;
    }

    char* substr(int pos) const {
        int length = size()-pos;
        return substr(pos, length);
    }

    void trim(int num_bases = 2){
        char *new_seq;
        int new_size = strlen(_seq) - 2*num_bases;

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

    void unassemble(){
        _assembled = false;
    }

    Read& operator = (const Read &rhs){
        _description = NULL;
        set_description(rhs._description);

        _seq = NULL;
        set_seq(rhs._seq);

        _qual = NULL;
        set_qual(rhs._qual);

        return *this;
    }

    ~Read(){
        if(_description != NULL){
            free(_description);
        }
        if(_seq != NULL){
            free(_seq);
        }
        if(_qual != NULL){
            free(_qual);
        }
    }
};

#endif
