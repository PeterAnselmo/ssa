#include <iostream>
#include <algorithm>
#include "assembly.cpp"
#include "fasta.cpp"

using namespace std;

string fastapath = "out.fasta";

int main(int argc, char* argv[]){

    if(argc != 2){
        printf("Please pass exactly one argument - the path to the input fastq file.\n");
        exit(1);
    }
    FastqFile fastq(argv[1]);
    printf("Fastq file Read.\n");
    if(DEBUGGING3){
        fastq.print_contents();
    }

    Assembly assem(fastq);
    assem.assemble_perfect_contigs();
    assem.trim_contigs();
    printf("Perfect Contigs assembled.\n");
    assem.print_contigs();
    assem.assemble_contigs();
    printf("Assembly completed, resulting in %u contigs.\n", static_cast<unsigned int>(assem.contigs.size()));
    assem.print_report();

    Fasta fasta(fastapath);
    fasta.description("SSA output");
    fasta.seq(assem.final_seq());
    fasta.write();
    printf("Fasta file Written.\n");

    return 0;
}

