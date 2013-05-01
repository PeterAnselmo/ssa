#ifndef FASTQ_CPP
#define FASTQ_CPP

#include <fstream>
#include <vector>
#include <stdlib.h>
#include "read.cpp"

using namespace std;

class Fastq {
private:
    vector<Read> _reads;

public:
    Fastq(char* filename){
        import_from_file(filename);
    }

    void import_from_file(char* filename){
        FILE *fp = fopen(filename, "r");
        if(fp == NULL){
            printf("Error Reading input Fastq file");
            exit(1);
        }
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

                //if read was worth keeping
                if(read.trim()){
                    _reads.push_back(read);
                }
            }
            ++count;
        }

        fclose(fp);
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
