#include <iostream>
#include <cassert>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samfile.cpp"

using namespace std;

int main(int argc, char* argv[]){
    string seq1 = "ACACACTATTGG";
    string seq2 = "ATGACANNNNNNAGCACACATTGG";

    cout << "Score: " << Assembly::smith_waterman(seq1, seq2) << endl;

    
    seq2 = "ACACACTATTGG";
    seq1 = "ATGACANNNNNNAGCACACATTGG";

    cout << "Score: " << Assembly::smith_waterman(seq1, seq2) << endl;


    seq1 = "ACACACTA";
    seq2 = "AGCACACA";

    assert(Assembly::smith_waterman(seq1, seq2) == 12);

}

