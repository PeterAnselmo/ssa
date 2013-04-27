#ifndef SW_MATRIX_H
#define SW_MATRIX_H

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "settings.cpp"

using namespace std;


class SWMatrix {
private:
    int **h;
    //int **s;
    int height;
    int width;
    char* _seq1; //along top of matrix
    char* _seq2; //along left side of matrix
    int max_h, max_w; //coordinate of the highest score in matrix
    char* gapped_seq1;
    char* gapped_seq2;

public:
    SWMatrix(char* new_seq1, char* new_seq2){

        //coincidentally, the matrix will have the same dimmensions as the allocated 
        //strings, since there is an extra zero row and a null character respectively.
        width = strlen(new_seq1) + 1;
        height = strlen(new_seq2) + 1;
        _seq1 = (char*)malloc(width);
        _seq2 = (char*)malloc(height);
        strncpy(_seq1, new_seq1, width);
        strncpy(_seq2, new_seq2, height);
        
        gapped_seq1 = NULL;
        gapped_seq2 = NULL;
        
        max_h = max_w = 0;

        if(DEBUGGING2){
            printf("Constructing matrix of wxh: %dx%d\n", width, height);
        }

        h = new int*[height];
        for(int i=0; i<height; ++i){
            h[i] = new int[width];
        }
        for(int i=0; i<height; ++i){
            h[i][0] = 0;
        }
        for(int j=0; j<width; ++j){
            h[0][j] = 0;
        }

        /*
        s = new int*[height];
        for(int i=0; i<height; ++i){
            s[i] = new int[width];
        }
        */

        compute_matrix();
        //compute_sum_matrix();
    }
    char* get_gapped_seq1(){
        return gapped_seq1;
    }
    char* get_gapped_seq2(){
        return gapped_seq2;
    }

    void compute_matrix(){
        const int w_match = 2;
        const int w_mismatch = -4;

        int max_val = 0;
        for(int i=1; i<height; i++){
            for(int j=1; j<width; ++j){
                int w = (_seq1[j-1] == _seq2[i-1]) ? w_match : w_mismatch;
                h[i][j] = max(
                    max(0,
                        h[i-1][j-1] + w), //match, mismatch
                    max(h[i-1][j] + w_mismatch, //deletion
                        h[i][j-1] + w_mismatch //insertion
                    )
                    );

                if(h[i][j] > max_val){
                    max_val = h[i][j];
                    max_h = i;
                    max_w = j;
                }
            }
        }
    }

    void print_matrix(){
        int padding = 3;
        cout << setw(padding) << " " << setw(padding) << '-';
        for(int i=0; i<width; ++i){
            cout << setw(padding+1) << _seq1[i];
        }
        cout << endl;
        cout << setw(padding-1) << '-';
        for(int i=0; i<width; ++i){
            cout << setw(padding+1) << '0';
        }
        cout << endl;
        for(int i=1; i<height; ++i){
            cout << setw(padding) << _seq2[i-1];
            for(int j=0; j<width; ++j){
                cout << setw(padding) << h[i][j] << " ";
            }
            cout << endl;
        }
    }

    int score(){
        //this will be the highest value in the matrix
        return h[max_h][max_w];
    }

    /*
    void compute_sum_matrix(){
        //start computation one column in, last column/row will always be 
        //own prefix max

        s[height-1][width-1] = h[height-1][width-1];
        //initialize last column
        for(int i=height-2; i>=0; --i){
            s[i][width-1] = h[i][width-1] + s[i+1][width-1];
        }
        //initialize last row
        for(int j=width-2; j>=0; --j){
            s[height-1][j] = h[height-1][j] + s[height-1][j+1];
        }
        for(int i=height-2; i>=0; --i){
            for(int j=width-2; j>=0; --j){
                s[i][j] = h[i][j] + max({
                    s[i+1][j], //down, insertion
                    s[i][j+1], //right, deletion
                    s[i+1][j+1], //diag, match/mismatch
                });
            }
        }
    }
    */

    /*
    void gap_seqs(){
        //these will hold the sequences with "-" in gaps
        //initially they will be reversed - more efficient to append,
        //then we will reverse them at the end
        
        //x & y will trace the path from lower right to top left
        int x = max_w;
        int y = max_h;
        while(x > 0 && y > 0){
            if(_seq1[x] == _seq2[y]){
                gapped_seq1 += _seq1[x];
                gapped_seq2 += _seq1[x];

            //if we're moving diagonal - match/mismatch
            if( h[y-1][x-1] >= h[y-1][x] && h[y-1][x-1] >= h[y][x-1]){
                //if(DEBUGGING3){ cout << 'a '; }
                gapped_seq1 += _seq1[x];
                gapped_seq2 += _seq2[y];
                --x;
                --y;

            //moving left - deletion in second sequence
            } else if ( h[y][x-1] > h[y-1][x] ){
                //if(DEBUGGING3){ cout << 'b ';}
                gapped_seq1 += _seq1[x];
                gapped_seq2 += '-';
                --x;

            //moving up - insertion in second sequence
            } else {
                //if(DEBUGGING3){ cout << 'c '; }
                gapped_seq1 += '-';
                gapped_seq2 += _seq2[y];
                --y;
            }
        }
        reverse(gapped_seq1.begin(), gapped_seq1.end());
        reverse(gapped_seq2.begin(), gapped_seq2.end());
    }
    */

    ~SWMatrix(){
        free(_seq1);
        free(_seq2);

        for(int i=0; i<height; ++i){
            delete[] h[i];
        }
        delete[] h;
    }

};

#endif
