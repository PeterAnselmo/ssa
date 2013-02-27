#include <iostream>
#include "fastafile.cpp"
#include "assembly.cpp"
#include "samfile.cpp"

using namespace std;

char outpath[] = "asdf.txt";

int main(int argc, char* argv[]){
    list<string> reads;

    FastaFile fasta(argv[1]);

    Assembly assem(fasta);
    assem.assemble();

    SamFile sam(assem);
    sam.write(outpath);

    cout << "hello world" << endl;

    return 0;
}

