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

    //get reds from input file
    Fastq input_file(argv[1]);
    printf("Fastq file Read.\n");
    if(DEBUGGING3){
        input_file.print_contents();
    }
    vector<Read> reads = input_file.reads();

    //Create assembly object, this will handle all the heavy lifting.
    Assembly assem(reads);

    //phase 1 - assemble perfect contigs
    assem.assemble_perfect_contigs();
    assem.trim_contigs();
    printf("Perfect Contigs assembled.\n");
    if(DEBUGGING){
        assem.print_contigs();
    }

    //phase 2 - merge contigs, allowing mismatches
    assem.assemble_contigs();
    printf("Assembly completed, resulting in %u contigs.\n", static_cast<unsigned int>(assem.contigs.size()));
    assem.print_report();

    //generate fasta file output
    if(assem.contigs.size() >= 1){
        Fasta fasta(fastapath);
        fasta.description("SSA output");
        fasta.seq(assem.final_seq());
        fasta.write();
        printf("Fasta file Written.\n");
    }

    return 0;
}

