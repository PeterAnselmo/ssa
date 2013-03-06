#ifndef FASTQFILE_H
#define FASTQFILE_H

#include <fstream>
#include <list>
#include <string>

using namespace std;

struct FQRead {
    string description;
    string seq;
    string plus;
    string qual;
};

class FastqFile {
public:
    list<FQRead> reads;

public:
    FastqFile(char* filename){
        import_from_file(filename);
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
        FQRead read;
        while(getline(fh, line)){
            if( count % 4 == 0 ){
                read.description = line;
            } else if (count % 4 == 1 ){
                read.seq = line;
            } else if (count % 4 == 2 ){
                read.plus = line;
            } else if (count % 4 == 3 ){
                read.qual = line;
                reads.push_back(read);
            }
            ++count;
        }

        fh.close();
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
