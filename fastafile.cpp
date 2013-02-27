#ifndef FASTAFILE_H
#define FASTAFILE_H

#include <fstream>
#include <list>
#include <string>

using namespace std;

class FastaFile {
public:
    list<string> reads;

public:
    FastaFile(char* filename){
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
        while(getline(fh, line)){
            reads.push_back(line);
        }

        fh.close();
    }
};

#endif
