#include <iostream>
#include <algorithm>
#include "assembly.cu"

using namespace std;

char outpath[] = "asdf.txt";
//char inpath[] = "/home/audioman/Storage/BioInfo/reads/glados_758/phix_10k.fastq";
char inpath[] = "/home/audioman/Storage/BioInfo/reads/illumina_phix/phix_10_by_hand.fastq";
string fastapath = "out.fasta";

int main(int argc, char* argv[]){

    FastqFile fastq(inpath);
    /*
    cout << "File Contents: " << endl;
    fastq.print_contents();
    */

    Assembly assem(fastq);
    assem.assemble_perfect_contigs();
//    assem.trim_contigs();
    assem.assemble_contigs();

    //sort the reads before displaying below contigs
    /*
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
        */

    Fasta fasta(fastapath);
    fasta.description("SSA output");
    fasta.seq(assepm.contigs.front().seq());
    fasta.write();

    return 0;
}

