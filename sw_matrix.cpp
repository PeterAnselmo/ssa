#ifndef SW_MATRIX_H
#define SW_MATRIX_H

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;

class SWMatrix {
private:
    int **h;
    int height;
    int width;
    string seq1; //along top of matrix
    string seq2; //along left side of matrix
    string gapped_seq1;
    string gapped_seq2;

public:
    SWMatrix(string new_seq1, string new_seq2){

        seq1 = new_seq1;
        seq2 = new_seq2;
        width = seq1.size() + 1;
        height = seq2.size() + 1;

        h = new int*[height];
        for(int i=0; i<height; ++i){
            h[i] = new int[width];
        }

        compute_matrix();
    }
    string get_gapped_seq1(){
        return gapped_seq1;
    }
    string get_gapped_seq2(){
        return gapped_seq2;
    }

    void compute_matrix(){
        const int w_match = 2;
        const int w_mismatch = -4;

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
    }

    void print_matrix(){
        int padding = 4;
        cout << setw(padding) << " " << setw(padding) << '-';
        for(int i=0; i<width; ++i){
            cout << setw(padding+1) << seq1[i];
        }
        cout << endl;
        cout << setw(padding-1) << '-';
        for(int i=0; i<width; ++i){
            cout << setw(padding+1) << '0';
        }
        cout << endl;
        for(int i=1; i<height; ++i){
            cout << setw(padding) << seq2[i-1];
            for(int j=0; j<width; ++j){
                cout << setw(padding) << h[i][j] << " ";
            }
            cout << endl;
        }
    }

    int score(){
        return h[height-1][width-1];
    }

    void gap_seqs(){
        //these will hold the sequences with "-" in gaps
        //initially they will be reversed - more efficient to append,
        //then we will reverse them at the end
        
        //x & y will trace the path from lower right to top left
        int x = width-1;
        int y = height-1;
        while(x > 0 && y > 0){
            /*
            if(seq1[x] == seq2[y]){
                gapped_seq1 += seq1[x];
                gapped_seq2 += seq1[x];
                */

            //if we're moving diagonal - match/mismatch
            if( h[y-1][x-1] >= h[y-1][x] && h[y-1][x-1] >= h[y][x-1]){
                cout << 'a' << endl;
                gapped_seq1 += seq1[x];
                gapped_seq2 += seq2[y];
                --x;
                --y;

            //moving left - deletion in second sequence
            } else if ( h[y][x-1] > h[y-1][x] ){

                cout << 'b' << endl;
                gapped_seq1 += seq1[x];
                gapped_seq2 += '-';
                --x;

            //moving up - insertion in second sequence
            } else {
                cout << 'c' << endl;
                gapped_seq1 += '-';
                gapped_seq2 += seq2[y];
                --y;
            }
        }
        reverse(gapped_seq1.begin(), gapped_seq1.end());
        reverse(gapped_seq2.begin(), gapped_seq2.end());
    }

    ~SWMatrix(){

        for(int i=0; i<height; ++i){
            delete[] h[i];
        }
        delete[] h;

    }

};

#endif
