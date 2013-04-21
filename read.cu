#ifndef READ_CU
#define READ_CU

#include <algorithm>
#include <string>

using namespace std;

class Read {
private:
    string _description;
    string _seq;
    string _plus;
    string _qual;
    bool _assembled;
public:
    unsigned int assem_contig;
    string gapped_seq;
    int assem_pos;

    __host__ __device__ Read(){
        _assembled = false;
    }

    string::size_type size() const {
        return _seq.size();
    }

    bool assembled() const{
        return _assembled;
    }
    void assemble(unsigned int contig_id, unsigned int pos){
        _assembled = true;
        assem_contig = contig_id;
        assem_pos = pos;
        gapped_seq = _seq;
    }
    void unassemble(){
        _assembled = false;
    }

    const string description() const{
        return _description;
    }
    void description(string new_description){
        _description = new_description;
    }

    const string seq() const{
        return _seq;
    }
    void seq(string new_seq){
        _seq = new_seq;
    }

    const string plus() const{
        return _plus;
    }
    void plus(string new_plus){
        _plus = new_plus;
    }

    const string qual() const{
        return _qual;
    }
    void qual(string new_qual){
        _qual = new_qual;
    }

    const string substr(int pos, int length) const {
        return _seq.substr(pos, length);
    }

    string rev_comp(){
        string rev_comp = "";

        for(int i=_seq.length()-1; i >= 0; --i){
            switch(_seq[i]){
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
        _seq = rev_comp();
        reverse(_qual.begin(), _qual.end());
    }

};

#endif
