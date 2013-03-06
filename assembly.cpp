#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <algorithm>
#include <list>
#include "fastqfile.cpp"

using namespace std;


class Assembly {
public:
    //someconatiner kmer graph

    //Established reference, used for glocal anchoring
    //string given_reference;
    
    //assembled reference
    //string reference;

    //std::list<string> unmatched_reads;

public:
    Assembly(FastqFile input_file){
    }
    void assemble(){

    }

    static int smith_waterman(string seq1, string seq2){
        const int w_match = 2;
        const int w_mismatch = -1;

        int width = seq1.size() + 1;
        int height = seq2.size() + 1;

        int h[height][width];

        for(int i=0; i<height; ++i){
            h[0][i] = 0;
        }
        for(int j=0; j<width; ++j){
            h[j][0] = 0;
        }
        for(int i=1; i<height; i++){
            for(int j=1; j<width; ++j){
                int w = (seq1[j-1] == seq2[i-1]) ? w_match : w_mismatch;
                h[i][j] = max({
                    0,
                    h[i-1][j-1] + w, //match, mismatch
                    h[i-1][j] + w_mismatch, //deletion
                    h[i][j-1] + w_mismatch //insertion
                    });
            }
        }
        /*
        for(int i=0; i<height; ++i){
            for(int j=0; j<width; ++j){
                cout << h[i][j] << " ";
            }
            cout << endl;
        }
        */
        return h[height-1][width-1];
    }

};

#endif
