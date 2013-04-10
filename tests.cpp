#include <iostream>
#include <cassert>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samfile.cpp"

using namespace std;

int main(int argc, char* argv[]){

    /*
    SWMatrix sw("ACACACTATTGG", "ATGACANNNNNNAGCACACATTGG");
    assert(sw.score() == 18);

    SWMatrix sw2("ATGACANNNNNNAGCACACATTGG", "ACACACTATTGG");
    assert(sw2.score() == 18);
    */

    SWMatrix sw3("ACACACTA", "AGCACACA");
    sw3.print_matrix();
    sw3.gap_seqs();
    cout << "Gapped Seqs:\n";
    cout << sw3.get_gapped_seq1() << endl;
    cout << sw3.get_gapped_seq2() << endl;
    assert(sw3.score() == 10);

    Read read;
    read.seq = "ATACGA";
    assert(read.rev_comp() == "TCGTAT");

    read.set_rev_comp();
    assert(read.seq == "TCGTAT");

    cout << "All tests passed." << endl;

}

