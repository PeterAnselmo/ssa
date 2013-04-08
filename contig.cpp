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
    string get_seq(){
        return seq;
    }
    void set_seq(string new_seq){
        seq = new_seq;
        qual = "";
        for(unsigned int i=0; i< seq.size(); ++i){
           qual += "1";
        }
    }
    string get_qual(){
        return qual;
    }
    
    string::size_type size(){
        return seq.size();
    }
    string substr(int start, int length){
        return seq.substr(start, length);
    }
    
    void append(string new_seq){
        seq += new_seq;
        for(unsigned int i=0; i< new_seq.size(); ++i){
           qual += "1";
        }
    }

    void inc_qual(int pos){
        string scores = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

        char current_score = qual[pos];
        if( current_score == 'Z' ){
            return;
        }
        string::size_type current_pos = scores.find(current_score);
        if( current_pos == string::npos ){
            cout << "Error incrementing quality value";
            exit(1);
        }
        qual[pos] = scores[++current_pos];
    }
};

#endif
