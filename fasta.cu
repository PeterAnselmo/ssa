#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "read.cu"

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

    void write(int line_width = 80){
        ofstream fh;
        fh.open(_filename.c_str());
        if( fh.fail() ){
            cout << "Error opening Fasta file for writing.\n";
            exit(1);
        }
        fh << ">" << _description << endl;

        int num_lines = ceil(_seq.size()/80.0);
        for(int i=0; i<num_lines; ++i){
            fh << _seq.substr(i*line_width, line_width) << endl;
        }
    }

};

#endif
