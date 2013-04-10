#ifndef READ_H
#define READ_H

#include <algorithm>
#include <string>

using namespace std;

class Read {
public:
    string description;
    string seq;
    string plus;
    string qual;
    int assem_contig;
    string gapped_seq;
    int assem_pos;

    const string::size_type size() const {
        return seq.size();
    }

    const string substr(int pos, int length) const {
        return seq.substr(pos, length);
    }

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
