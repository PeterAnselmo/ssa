#include <iostream>
#include <algorithm>
#include "assembly.cpp"
#include "fasta.cpp"

using namespace std;

char outpath[] = "asdf.txt";
char inpath[] = "/home/audioman/Storage/BioInfo/reads/glados_758/phix_10k.fastq";
//char inpath[] = "/home/audioman/Storage/BioInfo/reads/illumina_phix/phix_10_by_hand.fastq";
string fastapath = "out.fasta";

int main(int argc, char* argv[]){

    FastqFile fastq(inpath);
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
    printf("Merged Contigs assembled.\n");
    assem.print_report();

    Fasta fasta(fastapath);
    fasta.description("SSA output");
    fasta.seq(assem.contigs.front().seq());
    fasta.write();
    printf("Fasta file Written.\n");

    return 0;
}

