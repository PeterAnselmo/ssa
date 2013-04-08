#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <algorithm>
#include <list>
#include "contig.cpp"
#include "seqread.cpp"
#include "fastqfile.cpp"
#include <iomanip>

using namespace std;

const bool DEBUGGING = true;
const bool DEBUGGING2 = false;
const bool DEBUGGING3 = false;
//during sw this many bases on the right of the read will not be considered
//this is to remove the bias against aligning the read past the consensus
//const int TAIL_SIZE = 5;
const int MIN_OVERLAP = 15;
const int CONTIG_CAP = 1000;
const int MATCH_THRESHOLD = 40;

class Assembly {
public:
    list<SeqRead> reads;
    list<Contig> contigs;
    string reference;

public:
    Assembly(FastqFile &input_file){
        reads = input_file.reads;
        if(reads.empty()){
            cout << "[ERROR] no reads to align." << endl;
            exit(1);
        }

        reference = "";
    }
    void assemble(){

        int pass = 0;
        while(pass < CONTIG_CAP){
            Contig c(pass++);

            //find the first unmapped read
            bool all_mapped = true;
            for(auto &read : reads ){
                if(read.assem_pos == -1){
                    read.assem_pos = 0;
                    read.assem_contig = pass;
                    c.set_seq(read.seq);
                    all_mapped = false;
                    break;
                }
            }
            if( all_mapped ){
                return;
            }

            if(DEBUGGING) {
                cout << "Starting new contig with sequence: " << c.get_seq() << endl;
            }

            bool mapped_read = true; //if any reads were mapped in the last iteration
            while(mapped_read){
                mapped_read = false;

                if(DEBUGGING) {
                    cout << "Restarting at the beginning of read list\n";
                }
                list<SeqRead>::iterator iter;
                for( iter = reads.begin(); iter != reads.end(); ++iter){
                    SeqRead read = *iter;

                    //if this read was already mapped
                    if( read.assem_pos != -1 ){
                        continue;
                    }

                    cout << "Considering read: " << read.seq << endl;

                    int high_score = 0;
                    int high_pos = 0;
                    list<SeqRead>::iterator high_iter;

                    for(unsigned int i=0; i<c.size(); ++i){
                        int compare_size = min({read.seq.size(), c.size()-i});

                        if(compare_size < MIN_OVERLAP){
                            break;
                        }
                        int score = smith_waterman_score(c.substr(i,compare_size), read.seq.substr(0,compare_size));
                        int rev_score = smith_waterman_score(c.substr(i,compare_size), read.rev_comp().substr(0,compare_size));

                        if(DEBUGGING2){
                            cout << "Score at pos: " << i << ":" << score << endl;
                            cout << "RC Score at pos: " << i << ":" << rev_score << endl;
                        }
                        if( score > MATCH_THRESHOLD && score > high_score ){
                            high_score = score;
                            high_iter = iter;
                            high_pos = i;
                        }
                        if( rev_score > MATCH_THRESHOLD && rev_score > high_score ){
                            high_score = rev_score;
                            high_iter = iter;
                            high_pos = i;
                            read.set_rev_comp();
                        }
                    }
                    //ignoring assemblies at position 0, likely means that right half of seq
                    //aligned to before start of reference
                    if( high_score != 0 && high_pos != 0 ){
                        assemble_read(c, *high_iter, high_pos);
                        mapped_read = true;
                    }
                }
            }
            contigs.push_back(c);
        }
    }

    static int smith_waterman_score(string seq1, string seq2){
        int width = seq1.size() + 1;
        int height = seq2.size() + 1;
        if(DEBUGGING3){
            cout << "Computing SW Score " << seq1 << ":" << seq2 << endl;
        }

        int **h;
        h = new int*[height];
        for(int i=0; i<height; ++i){
            h[i] = new int[width];
        }

        smith_waterman_matrix(seq1, seq2, h, height, width);

        int score = h[height-1][width-1];
        if(DEBUGGING3){
            cout << "Score: " << score << endl;
        }

        for(int i=0; i<height; ++i){
            delete[] h[i];
        }
        delete[] h;

        return score;
    }

    void print_report(){
        int num_assembled = 0;
        for(auto &read : reads){
            if( read.assem_pos != -1 ){
                ++num_assembled;
            }
        }
        cout << "Assembled " << num_assembled << "/" << reads.size() << " reads ("
             << static_cast<double>(num_assembled) * 100 /reads.size() << "%)" << endl;

    }

private:
    void assemble_read(Contig &c, SeqRead &read, int pos){
        if(DEBUGGING){
            cout << "Assembling Read: " << read.seq << " to contig " << c.id << " at " << pos << endl;
        }
        read.assem_pos = pos;
        read.assem_contig = c.id;

        //if the read extends to the right of the c.seq
        unsigned int overlap_size = c.size() - pos; 
        if(DEBUGGING){
            cout << "overlap size: " << overlap_size << endl;
            cout << "Read size: " << read.seq.size() << endl;
        }
        if( overlap_size < read.seq.size() ){
            for(unsigned int i=0; i<overlap_size; ++i){
                c.inc_qual(pos+i);
            }
            string new_seq = read.seq.substr(overlap_size);
            if(DEBUGGING){
                cout << "adding " << new_seq << " to reference.\n";
            }
            c.append(new_seq);
        }
    }

    static void smith_waterman_matrix(const string &seq1, const string &seq2, int **h, int height, int width){
        const int w_match = 2;
        const int w_mismatch = -2;

        for(int i=0; i<height; ++i){
            h[i][0] = 0;
        }
        for(int j=0; j<width; ++j){
            h[0][j] = 0;
        }
        for(int i=1; i<height; i++){
            for(int j=1; j<width; ++j){
                int w = (seq1[j-1] == seq2[i-1]) ? w_match : w_mismatch;
                h[i][j] = max({
                    0,
                    h[i-1][j-1] + w, //match, mismatch
                    h[i-1][j] + w_mismatch, //deletion
                    h[i][j-1] + w_mismatch //insertion
                    });
            }
        }
        if(DEBUGGING3){
            for(int i=0; i<height; ++i){
                for(int j=0; j<width; ++j){
                    cout << setw(5) << h[i][j] << " ";
                }
                cout << endl;
            }
        }
    }

};

#endif
