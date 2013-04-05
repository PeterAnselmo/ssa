#ifndef CONGTIG_H
#define CONGTIG_H

#include <algorithm>
#include <string>

using namespace std;

class Contig {
public:
    int id;
    string seq;
    string qual;

    Contig(){}
    Contig(int new_id){id = new_id;}
};

#endif
