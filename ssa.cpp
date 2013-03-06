#include <iostream>
#include "fastqfile.cpp"
#include "assembly.cpp"
#include "samfile.cpp"

using namespace std;

char outpath[] = "asdf.txt";

int main(int argc, char* argv[]){
    list<string> reads;

    FastqFile fastq(argv[1]);
    cout << "File Contents: " << endl;
    fastq.print_contents();

    Assembly assem(fastq);
    assem.assemble();

    SamFile sam(assem);
    sam.write(outpath);

    return 0;
}

