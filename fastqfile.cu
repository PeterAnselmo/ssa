#ifndef FASTQFILE_CU
#define FASTQFILE_CU

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "read.cu"

using namespace std;

//number of bases to delete off both ends of the read
const int TRIM_SIZE = 2;

class FastqFile {
private:
    Read *_reads;
    int num_reads;

public:
    FastqFile(char* filename){
        //todo, make this smarter and detect num of reads
        _reads = new Read[25000];
        import_from_file(filename);
        trim_reads();
    }

    void import_from_file(char* filename){
        ifstream fh;
        fh.open(filename);
        FILE *fp = fopen(filename, "r");
        char line[256];

        int count = 0;
        int read_num = 0;
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
                ++read_num;
                read.set_qual(line);
                Read *temp = new Read(read);
                _reads[read_num] = *temp;
            }
            ++count;
        }

        fh.close();
    }

    void trim_reads(){
        for(int i=0; i<num_reads; ++i){
            _reads[i].trim();
        }
    }

    Read* reads() const{
        return _reads;
    }

    void print_contents(){
        for(int i=0; i<num_reads; ++i){
            printf("%s\n%s\n%c\n%s\n", 
                    _reads[i].description(), 
                    _reads[i].seq(), 
                    '+',
                    _reads[i].qual());
        }
    }
};

#endif
