#ifndef FASTQFILE_CU
#define FASTQFILE_CU

#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "read.cpp"

using namespace std;

//number of bases to delete off both ends of the read
const int TRIM_SIZE = 2;

class FastqFile {
private:
    vector<Read> _reads;

public:
    FastqFile(char* filename){
        import_from_file(filename);
        trim_reads();
    }

    void import_from_file(char* filename){
        ifstream fh;
        fh.open(filename);
        FILE *fp = fopen(filename, "r");
        char line[256];

        int count = 0;
        Read read;
        while(fgets(line,sizeof(line), fp)){
            //chomp newline
            line[strlen(line)-1] = '\0';

            if( count % 4 == 0 ){
                read.set_description(line);
            } else if (count % 4 == 1 ){
                read.set_seq(line);
                read.set_gapped_seq(line);
            } else if (count % 4 == 2 ){
                //read.plus(line.c_str());
            } else if (count % 4 == 3 ){
                read.set_qual(line);
                _reads.push_back(read);
            }
            ++count;
        }

        fh.close();
    }

    void trim_reads(){
        vector<Read>::iterator read;
        for(read = _reads.begin(); read != _reads.end(); ++read){
            read->trim();
        }
    }

    const vector<Read> reads() const{
        return _reads;
    }

    void print_contents(){
        vector<Read>::iterator read;
        for(read = _reads.begin(); read != _reads.end(); ++read){
            printf("%s\n%s\n%c\n%s\n", 
                    read->description(), 
                    read->seq(), 
                    '+',
                    read->qual());
        }
    }
};

#endif
