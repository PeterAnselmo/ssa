#ifndef SAMFILE_H
#define SAMFILE_H
#include "assembly.cpp"
#include "samtools-0.1.19/sam.h"

using namespace std;

class SamFile {
private:
    string _filepath;
    
public:
    SamFile(string new_filepath, Assembly assem){
        _filepath = new_filepath;
        write();
    }
    void write(){

        samfile_t *fh = NULL;
        fh = samopen(_filepath.c_str(), "w", 0);

        if(NULL == fh){
           printf("Could not read sam file"); 
        }
        
        bam1_core_t core;
        core.tid=1;
        core.pos = 2;
        core.l_qseq = 3;

        bam1_t *bam;
        bam->core = core;

        int written = samwrite(fh, bam);
    }
};

#endif
