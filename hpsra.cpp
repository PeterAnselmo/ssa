#include <iostream>
#include <list>
#include <fstream>

using namespace std;

void populate_reads(list<string> &, char[]);

int main(int argc, char* argv[]){
    list<string> reads;

    populate_reads(reads, argv[1]);

    for(const auto &read : reads){
        cout << read << endl;
    }
    cout << "hello world" << endl;

    return 0;
}

void populate_reads(list<string> &reads, char* filename){

    ifstream fh;
    fh.open(filename);

    if( fh.fail() ){
        cout << "Error opening read file.\n";
        exit(1);
    }

    string line;
    while(getline(fh, line)){
        reads.push_back(line);
    }

    fh.close();

}
