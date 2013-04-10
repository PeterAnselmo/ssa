#ifndef CONGTIG_H
#define CONGTIG_H

#include <algorithm>
#include <string>

using namespace std;

class Contig {
private:
    string seq;
    string qual;
public:
    int id;

    Contig(){}
    Contig(int new_id){id = new_id;}
    Contig(int new_id, string new_seq){
        id = new_id;
        set_seq(new_seq);
    }

    const int get_id() const{
        return id;
    }

    const string get_seq() const{
        return seq;
    }

    void set_seq(string new_seq){
        seq = new_seq;
        qual = "";
        for(unsigned int i=0; i< seq.size(); ++i){
           qual += "!";
        }
    }

    const string get_qual() const {
        return qual;
    }
    
    const string::size_type size(){
        return seq.size();
    }
    const string substr(int start, int length){
        return seq.substr(start, length);
    }
    
    void append(string new_seq){
        seq += new_seq;
        for(unsigned int i=0; i< new_seq.size(); ++i){
           qual += "!";
        }
    }

    void inc_qual(int pos){
        qual[pos] = qual[pos] + 1;
    }
};

#endif
