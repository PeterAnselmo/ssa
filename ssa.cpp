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
    cout << "Fastq file Read." << endl;
    fastq.print_contents();

    Assembly assem(fastq);
    assem.assemble_perfect_contigs();
    assem.trim_contigs();
//    assem.assemble_contigs();
    assem.print_contigs();
    cout << "Contigs assembled." << endl;

    Fasta fasta(fastapath);
    fasta.description("SSA output");
    fasta.seq(assem.contigs.front().seq());
    fasta.write();
    cout << "Fasta file Written." << endl;

    return 0;
}

