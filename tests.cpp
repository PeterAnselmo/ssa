#include <iostream>
#include <cassert>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samfile.cpp"

using namespace std;

int main(int argc, char* argv[]){
    string seq1 = "ACACACTA";
    string seq2 = "AGCACACA";

    assert(Assembly::smith_waterman(seq1, seq2) == 12);
}

