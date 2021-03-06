#ifndef READ_CU
#define READ_CU

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include "settings.cpp"
#include "utils.cpp"

class Read {
private:
    
    char* _description;
    char* _seq;
    char* _rev_comp;
    char* _qual;
    char* _gapped_seq;

    bool _assembled;
    unsigned int _contig;
    unsigned int _position;
public:

    Read(){
        _description = NULL;
        _seq = NULL;
        _rev_comp = NULL;
        _qual = NULL;
        _gapped_seq = NULL;

        _assembled = false;
        _contig = 0;
        _position = 0;
    }
    Read(const Read &rhs){
        _description = NULL;
        _seq = NULL;
        _rev_comp = NULL;
        _qual = NULL;
        _gapped_seq = NULL;

        set_description(rhs._description);
        set_seq(rhs._seq, false);
        set_rev_comp(rhs._rev_comp);
        set_qual(rhs._qual);
        set_gapped_seq(rhs._gapped_seq);

        _assembled = rhs._assembled;
        _contig = rhs._contig;
        _position = rhs._position;
    }

    bool assembled() const{
        return _assembled;
    }

    void assemble(unsigned int contig_id, unsigned int pos){
        _assembled = true;
        _contig = contig_id;
        _position = pos;
        set_gapped_seq(_seq);
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

    char* gapped_seq(){
        return _gapped_seq;
    }

    unsigned int gapped_size() const {
        return strlen(_gapped_seq);
    }

    char* gapped_substr(int pos, int length) const {
        char *sub = (char*)malloc(length+1);
        check_ptr(sub);
        memcpy(sub, &_gapped_seq[pos], length);
        sub[length] = '\0';
        return sub;
    }
    char* gapped_substr(int pos) const {
        int length = strlen(_gapped_seq) - pos;
        return gapped_substr(pos, length);
    }

    void init_gapped_seq(){
        int length = strlen(_seq)+1;
        if(_gapped_seq == NULL){
            _gapped_seq = (char*)malloc(length);
        } else {
            _gapped_seq = (char*)realloc(_gapped_seq, length);
        }
        memcpy(_gapped_seq, &_seq[0], length);
    }

    unsigned int position(){
        return _position;
    }

    const char* qual() const{
        return _qual;
    }

    char* qual_substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        check_ptr(sub);

        memcpy(sub, &_qual[pos], length);
        sub[length] = '\0';
        return sub;
    }

    char* rev_substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        check_ptr(sub);

        memcpy(sub, &_rev_comp[pos], length);
        sub[length] = '\0';
        return sub;
    }

    char* rev_comp() const{
        return _rev_comp;
    }

    const char* seq() const{
        return _seq;
    }

    void compute_rev_comp(){
        if(_rev_comp == NULL){
            _rev_comp = (char*)malloc(size()+1);
        } else {
            _rev_comp = (char*)realloc(_rev_comp, size()+1);
        }
        check_ptr(_rev_comp);

        for(unsigned int i=0; i <size(); ++i){
            int pos = size()-1-i;
            switch(_seq[pos]){
                case 'A':
                    _rev_comp[i] = 'T';
                break;
                case 'T':
                    _rev_comp[i] = 'A';
                break;
                case 'C':
                    _rev_comp[i] = 'G';
                break;
                case 'G':
                    _rev_comp[i] = 'C';
                break;
                case 'N':
                    _rev_comp[i] = 'N';
                break;
            }
        }
        _rev_comp[size()] = '\0';
    }

    void set_description(const char *new_description){
        if(_description == NULL){
            _description = (char*)malloc(strlen(new_description)+1);
        } else {
            _description = (char*)realloc(_description, strlen(new_description)+1);
        }
        check_ptr(_description);

        strncpy(_description, new_description, strlen(new_description)+1);
    }

    void set_gapped_seq(const char* new_gapped_seq){
        if(_gapped_seq == NULL){
            _gapped_seq = (char*)malloc(strlen(new_gapped_seq)+1);
        } else {
            _gapped_seq = (char*)realloc(_gapped_seq, strlen(new_gapped_seq)+1);
        }
        check_ptr(_gapped_seq);

        strncpy(_gapped_seq, new_gapped_seq, strlen(new_gapped_seq)+1);
    }

    void set_position(unsigned int new_pos){
        _position = new_pos;
    }

    void set_qual(const char* new_qual){
        if(_qual == NULL){
            _qual = (char*)malloc(strlen(new_qual)+1);
        } else {
            _qual = (char*)realloc(_qual, strlen(new_qual)+1);
        }
        check_ptr(_qual);

        strncpy(_qual, new_qual, strlen(new_qual)+1);
    }

    void set_rev_comp(const char* new_rev_comp){
        if(_rev_comp == NULL){
            _rev_comp = (char*)malloc(strlen(new_rev_comp)+1);
        } else {
            _rev_comp = (char*)realloc(_rev_comp, strlen(new_rev_comp)+1);
        }
        check_ptr(_rev_comp);

        strncpy(_rev_comp, new_rev_comp, strlen(new_rev_comp)+1);
    }

    void set_seq(const char* new_seq, bool set_rev_comp = true){
        int length = strlen(new_seq);
        if(_seq == NULL){
            _seq = (char*)malloc(length+1);
        } else {
            _seq = (char*)realloc(_seq, length+1);
        }
        check_ptr(_seq);

        strncpy(_seq, new_seq, length+1);
        //_seq[length] = '\0';
        if( set_rev_comp ){
            compute_rev_comp();
        }
    }

    unsigned int size() const {
        return strlen(_seq);
    }

    char* substr(int pos, int length) const {
        char *sub = (char *)malloc(length+1);
        check_ptr(sub);

        memcpy(sub, &_seq[pos], length);
        sub[length] = '\0';
        return sub;
    }

    char* substr(int pos) const {
        int length = size()-pos;
        return substr(pos, length);
    }

    void swap_rev_comp(){
        char* temp = _seq;
        _seq = _rev_comp;
        _rev_comp = temp;
    }

    bool trim(int num_bases = 2){
        char *new_seq;
        char *new_qual;

        //trim left side
        unsigned int start_pos = num_bases; 
        while(start_pos < size() && _seq[start_pos] == 'N'){
            ++start_pos;
        }

        //trim right side
        unsigned int end_pos = size()-num_bases;
        while(end_pos > 0 && _seq[end_pos-1] == 'N'){
            --end_pos;
        }

        int new_size = end_pos - start_pos;
        //if read is too short to use after trimming desired amount
        if( new_size <= MIN_OVERLAP){
            return false;
        }


        new_seq = substr(start_pos, new_size);
        new_qual = qual_substr(start_pos, new_size);
        free(_seq);
        free(_qual);
        _seq = new_seq;
        _qual = new_qual;
        compute_rev_comp();
        return true;
    }

    void unassemble(){
        _assembled = false;
    }

    const bool operator < (const Read &rhs) const {

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
        set_seq(rhs._seq, false);
        set_rev_comp(rhs._rev_comp);
        set_qual(rhs._qual);
        set_gapped_seq(rhs._gapped_seq);

        _assembled = rhs._assembled;
        _contig = rhs._contig;
        _position = rhs._position;

        return *this;
    }

    ~Read(){
        free_strings();
    }

private:
    void free_strings(){
        if(_description != NULL){
            free(_description);
            _description = NULL;
        }
        if(_seq != NULL){
            free(_seq);
            _seq = NULL;
        }
        if(_rev_comp != NULL){
            free(_rev_comp);
            _rev_comp = NULL;
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
