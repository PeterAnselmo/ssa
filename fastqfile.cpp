#ifndef FASTQFILE_H
#define FASTQFILE_H

#include <fstream>
#include <list>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "read.cpp"

using namespace std;

//number of bases to delete off both ends of the read
const int TRIM_SIZE = 2;

class FastqFile {
private:
    list<Read> _reads;

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
        Read read;
        while(getline(fh, line)){
            if( count % 4 == 0 ){
                read.description(line);
            } else if (count % 4 == 1 ){
                read.seq(line);
            } else if (count % 4 == 2 ){
                read.plus(line);
            } else if (count % 4 == 3 ){
                read.qual(line);
                read.assem_pos = -1;
                _reads.push_back(read);
            }
            ++count;
        }

        fh.close();
    }

    void trim_reads(){
        for(auto &read : _reads ){
            read.seq(read.seq().substr(TRIM_SIZE, read.seq().size()-2*TRIM_SIZE));
        }
    }

    const list<Read> reads() const{
        return _reads;
    }

    void print_contents(){
        for(const auto read : _reads ){
            printf("%s\n%s\n%s\n%s\n", read.description().c_str(), 
                    read.seq().c_str(), 
                    read.plus().c_str(), 
                    read.qual().c_str());
        }
    }
};

#endif
