#ifndef FASTQFILE_H
#define FASTQFILE_H

#include <fstream>
#include <list>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "seqread.cpp"

using namespace std;

class FastqFile {
public:
    list<SeqRead> reads;

public:
    FastqFile(char* filename){
        import_from_file(filename);
        trim_reads();
    }
    void import_from_file(char* filename){
        ifstream fh;
        fh.open(filename);

        if( fh.fail() ){
            cout << "Error opening read file.\n";
            exit(1);
        }

        string line;
        int count = 0;
        SeqRead read;
        while(getline(fh, line)){
            if( count % 4 == 0 ){
                read.description = line;
            } else if (count % 4 == 1 ){
                read.seq = line;
            } else if (count % 4 == 2 ){
                read.plus = line;
            } else if (count % 4 == 3 ){
                read.qual = line;
                read.assembled_pos = -1;
                reads.push_back(read);
            }
            ++count;
        }

        fh.close();
    }

    //for now, just take off first two bases, and only read 50
    void trim_reads(){
        for(auto &elem : reads ){
            elem.seq = elem.seq.substr(2,30);
        }
    }

    void print_contents(){
        for(const auto read : reads ){
            printf("%s\n%s\n%s\n%s\n", read.description.c_str(), 
                    read.seq.c_str(), 
                    read.plus.c_str(), 
                    read.qual.c_str());
        }
    }
};

#endif
