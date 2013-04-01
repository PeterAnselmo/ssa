#include <iostream>
#include <cassert>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samfile.cpp"

using namespace std;

int main(int argc, char* argv[]){
    string seq1 = "ACACACTATTGG";
    string seq2 = "ATGACANNNNNNAGCACACATTGG";
    assert(Assembly::smith_waterman(seq1, seq2) == 20);

    seq2 = "ACACACTATTGG";
    seq1 = "ATGACANNNNNNAGCACACATTGG";
    assert(Assembly::smith_waterman(seq1, seq2) == 20);

    seq1 = "ACACACTA";
    seq2 = "AGCACACA";

    assert(Assembly::smith_waterman(seq1, seq2) == 12);

    SeqRead read;
    read.seq = "ATACGA";
    assert(read.rev_comp() == "TCGTAT");

    read.set_rev_comp();
    assert(read.seq == "TCGTAT");

    cout << "All tests passed." << endl;

}

