#ifndef SEQREAD_H
#define SEQREAD_H

#include <algorithm>
#include <string>

using namespace std;

class SeqRead {
public:
    string description;
    string seq;
    string plus;
    string qual;
    int assem_contig;
    int assem_pos;

    string rev_comp(){
        string rev_comp = "";

        for(int i=seq.length()-1; i >= 0; --i){
            switch(seq[i]){
                case 'A':
                    rev_comp += "T";
                break;
                case 'T':
                    rev_comp += "A";
                break;
                case 'C':
                    rev_comp += "G";
                break;
                case 'G':
                    rev_comp += "C";
                break;
                case 'N':
                    rev_comp += "N";
                break;
            }
        }
        return rev_comp;
    }

    void set_rev_comp(){
        seq = rev_comp();
        reverse(qual.begin(), qual.end());
    }

};

#endif
