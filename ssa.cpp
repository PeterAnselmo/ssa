#include <iostream>
#include <algorithm>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samtools-0.1.19/sam.h"

using namespace std;

char outpath[] = "asdf.txt";
char inpath[] = "/home/audioman/Storage/BioInfo/reads/glados_758/phix_100k.fastq";
//char inpath[] = "/home/audioman/Storage/BioInfo/reads/illumina_phix/phix_10_by_hand.fastq";

int main(int argc, char* argv[]){

    FastqFile fastq(inpath);
    /*
    cout << "File Contents: " << endl;
    fastq.print_contents();
    */

    Assembly assem(fastq);
    assem.assemble_perfect_contigs();
    assem.trim_contigs();
//    assem.assemble_contigs();

    //sort the reads before displaying below contigs
    assem.reads.sort([](const Read &r1, const Read &r2){ return r1.assem_pos < r2.assem_pos; });

    for(const auto &contig : assem.contigs){
        cout << "Assembled Contig " << contig.id() << ":\n" << contig.seq() << endl;
        cout << contig.qual() << endl;
        for(const auto &read : assem.reads){
            if(read.assembled() && read.assem_contig == contig.id()){
                for(int i=0; i<read.assem_pos; ++i){
                    cout << " ";
                }
                cout << read.gapped_seq << endl;
            }
        }
    }

    assem.print_report();
    cout << endl << endl;
    
    /*
    for(const auto &read : assem.reads){
        printf("read- seq:%s | c: %d | pos: %d\n", read.gapped_seq.c_str(), read.assem_contig, read.assem_pos);
    }
    */
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

