#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <list>
#include "fastafile.cpp"

using namespace std;

const int OPEN_GAP_PENALTY = 10;
const int EXTEND_GAP_PENALTY = 1;

class Assembly {
public:
    //someconatiner kmer graph

    //Established reference, used for glocal anchoring
    string given_reference;
    
    //assembled reference
    string reference;

    std::list<string> unmatched_reads;

public:
    Assembly(FastaFile input){
        unmatched_reads = input.reads;
    }
    void assemble(){

    }
    

};

#endif
