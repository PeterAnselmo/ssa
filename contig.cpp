#ifndef CONGTIG_H
#define CONGTIG_H

#include <algorithm>
#include <string>
#include <vector>
#include "read.cpp"

using namespace std;

class Contig {
private:
    string _seq;
    string _qual;
    unsigned int _id;

public:

    Contig(){}
    Contig(int id){_id = id;}
    Contig(int id, string new_seq){
        _id = id;
        seq(new_seq);
    }

    unsigned int id() const{
        return _id;
    }

    string seq() const{
        return _seq;
    }

    void seq(string seq){
        _seq = seq;
        _qual = "";
        for(unsigned int i=0; i< _seq.size(); ++i){
           _qual += "!";
        }
    }

    string qual() const {
        return _qual;
    }
    
    unsigned int size(){
        return _seq.size();
    }

    string substr(int start, int length){
        return _seq.substr(start, length);
    }
    
    void prepend(string seq){
        _seq = seq + _seq;
        for(unsigned int i=0; i< seq.size(); ++i){
           _qual = "!" + _qual;
        }
    }

    void append(string seq){
        _seq += seq;
        for(unsigned int i=0; i< seq.size(); ++i){
           _qual += "!";
        }
    }

    void trim(char min_quality, vector<Read> &reads){
        min_quality += 33;

        //trim_left_size
        unsigned int pos = 0;
        while(pos < _seq.size() && _qual[pos] < min_quality){
            ++pos;
        }
        _seq = _seq.substr(pos);
        _qual = _qual.substr(pos);
        shift_aligned_reads(pos, reads);

        //trim_right_size
        pos = _seq.size();
        while(pos > 0 && _qual[pos-1] < min_quality){
            --pos;
        }
        _seq = _seq.substr(0,pos);
        _qual = _qual.substr(0,pos);

    }

    void inc_qual(int pos){
        if(_qual[pos] < 126) {
            _qual[pos] = _qual[pos] + 1;
        }
    }

    void unshift_aligned_reads(unsigned int distance, vector<Read> &reads){
        vector<Read>::iterator read;
        for(read = reads.begin(); read != reads.end(); ++read){
            if(read->assembled() && read->contig() == _id){
                read->set_position(read->position() + distance);
            }
        }
    }
    void shift_aligned_reads(unsigned int distance, vector<Read> &reads){
        vector<Read>::iterator read;
        for(read = reads.begin(); read != reads.end(); ++read){
            if(read->assembled() && read->contig() == _id){
                int overlap = distance - read->position();
                if( overlap > 0){
                    char *gapped_substr = read->gapped_substr(overlap);
                    strcpy(read->gapped_seq, gapped_substr);
                    free(gapped_substr);
                } else {
                    read->set_position(read->position() - distance);
                }
            }
        }
    }

    void merge(Contig other){
        

    }
};

#endif