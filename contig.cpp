#ifndef CONGTIG_H
#define CONGTIG_H

#include <algorithm>
#include <string>

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

    const unsigned int id() const{
        return _id;
    }

    const string seq() const{
        return _seq;
    }

    void seq(string seq){
        _seq = seq;
        _qual = "";
        for(unsigned int i=0; i< _seq.size(); ++i){
           _qual += "!";
        }
    }

    const string qual() const {
        return _qual;
    }
    
    const string::size_type size(){
        return _seq.size();
    }

    const string substr(int start, int length){
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

    void inc_qual(int pos){
        _qual[pos] = _qual[pos] + 1;
    }

    void unshift_aligned_reads(unsigned int distance, list<Read> &reads){
        for(auto &read : reads){
            if(read.assembled() && read.assem_contig == _id){
                read.assem_pos += distance;
            }
        }
    }
};

#endif
