#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include "settings.cpp"

using namespace std;

class Fasta{
private:
    string _filename;
    string _description;
    string _seq;

public:
    Fasta(string new_filename){
        _filename = new_filename;
    }

    void description(string new_description){
        _description = new_description;
    }
    void seq(string new_seq){
        _seq = new_seq;
    }

    void write(){
        FILE *fh = fopen(_filename.c_str(), "w");
        if( fh == NULL ){
            printf("Error opening Fasta file for writing.\n");
            exit(1);
        }
        fprintf(fh, ">%s\n", _description.c_str());

        int num_lines = ceil(static_cast<double>(_seq.size())/FASTA_LINE_WIDTH);

        for(int i=0; i<num_lines; ++i){
            fprintf(fh, "%s\n", _seq.substr(i * FASTA_LINE_WIDTH, FASTA_LINE_WIDTH).c_str());
        }
    }

};

#endif
