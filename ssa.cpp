#include <iostream>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samtools-0.1.19/sam.h"

using namespace std;

char outpath[] = "asdf.txt";
char inpath[] = "/home/audioman/Storage/BioInfo/reads/phix_100k.fastq";

int main(int argc, char* argv[]){
    list<string> reads;

    FastqFile fastq(inpath);
    /*
    cout << "File Contents: " << endl;
    fastq.print_contents();
    */

    Assembly assem(fastq);
    assem.assemble();

    for(auto &contig : assem.contigs){
        cout << "Assembled Contig " << contig.id << ":\n" << contig.get_seq() << endl;
        cout << contig.get_qual() << endl;
    }

    assem.print_report();
    //cout << "\n\nReference:\n" << assem.reference;

    /*
    samfile_t *fp_in = NULL;
    fp_in = samopen(inpath, "r", 0);

    if(NULL == fp_in){
       printf("Could not read sam file"); 
    }
    */

    return 0;
}

