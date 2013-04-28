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
    int height;
    int width;
    char* _seq1; //along top of matrix
    char* _seq2; //along left side of matrix
    int max_h, max_w; //coordinate of the highest score in matrix
    char* _gapped_seq1;
    char* _gapped_seq2;
    char* _merged_seq;

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
        
        _gapped_seq1 = NULL;
        _gapped_seq2 = NULL;
        _merged_seq = NULL;
        
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

        compute_matrix();
        if(DEBUGGING3){
            print_matrix();
        }
    }
    char* get_gapped_seq1(){
        return _gapped_seq1;
    }
    char* get_gapped_seq2(){
        return _gapped_seq2;
    }
    char* merged_seq(){
        return _merged_seq;
    }

    void compute_matrix(){
        max_h = max_w = 0;

        int max_val = 0;
        for(int i=1; i<height; i++){
            for(int j=1; j<width; ++j){
                int w = (_seq1[j-1] == _seq2[i-1]) ? SW_W_MATCH : SW_W_MISMATCH;
                h[i][j] = max(
                    max(0,
                        h[i-1][j-1] + w), //match, mismatch
                    max(h[i-1][j] + SW_W_MISMATCH, //deletion
                        h[i][j-1] + SW_W_MISMATCH //insertion
                    )
                    );

                //if were' on the last row or column, look for a maximum
                if( (i==height-1 || j==width-1) && (h[i][j] > max_val)){
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

    //Traceback the path from the bottom right to the top left of the matrix
    void gap_seqs(){
        //these will hold the sequences with "-" in gaps
        //we'll allocate the worst case scenario (length of both strings) to avoid the need to reallocate
        _gapped_seq1 = (char*)malloc(width+height+1);
        _gapped_seq2 = (char*)malloc(width+height+1);
        _merged_seq = (char*)malloc(width+height+1);

        //this will be the position of the gapped sequence as we assign
        int current_position = width+height;
        
        //x & y will trace the path from lower right to top left
        //currently gives preferences as follows: match/mismatch, deletion, insertion
        //an improvment would be to use the one with the highest prev score
        int x = max_w;
        int y = max_h;

        if( DEBUGGING2){
            printf("Computing Gapped Sequences, starting with %d,%d...\n", y,x);
        }

        while(x > 0 && y > 0){

            //diagonal - match or mismatch
            //note, string position is -1 from grid position
            if( ((_seq1[x-1] == _seq2[y-1]) /*&& (h[y][x] == (h[y-1][x-1] + SW_W_MATCH))*/) ||
                ((_seq1[x-1] != _seq2[y-1]) && (h[y][x] == (h[y-1][x-1] + SW_W_MISMATCH))) ) {

                if(DEBUGGING3){ printf("Traceback: Diagonal\n"); }

                _gapped_seq1[current_position] = _seq1[x];
                _gapped_seq2[current_position] = _seq2[y];
                --x;
                --y;

            //moving left - deletion in second sequence
            } else if ( h[y][x] == (h[y][x-1] + SW_W_MISMATCH) ){
                
                if(DEBUGGING3){ printf("Traceback: Left\n");}

                _gapped_seq1[current_position] = _seq1[x];
                _gapped_seq2[current_position] = '-';
                --x;

            //moving up - insertion in second sequence
            } else if ( h[y][x] == (h[y-1][x] + SW_W_MISMATCH) ) {

                if(DEBUGGING3){ printf("Traceback: Up\n"); }

                _gapped_seq1[current_position] = '-';
                _gapped_seq2[current_position] = _seq2[y];
                --y;
            } else if( h[y][x] == 0 ){
                
                if(DEBUGGING3){ printf("Traceback Lost Path: Guessing Diagonal\n"); }

                _gapped_seq1[current_position] = _seq1[x];
                _gapped_seq2[current_position] = _seq2[y];
                --x;
                --y;
                
            } else {
                printf("Error, could not compute valid traceback");
                exit(1);
            }
            --current_position;
        }
        //Fencepost
        _gapped_seq1[current_position] = _seq1[x];
        _gapped_seq2[current_position] = _seq2[y];

        int new_length = width + height - current_position;
        memmove(&_gapped_seq1[0], &_gapped_seq1[current_position], new_length);
        _gapped_seq1[new_length] = '\0';
        memmove(&_gapped_seq2[0], &_gapped_seq2[current_position], new_length);
        _gapped_seq2[new_length] = '\0';

        if(DEBUGGING3){
            printf("Gapped Contigs:\n%s\n%s\n", _gapped_seq1, _gapped_seq2);
        }

        //1. assemble merge before match
        int pre_length;
        //seq1 extended to left of seq2
        if(x > 0 && y == 0){
            pre_length = x;
            memcpy(&_merged_seq[0], &_seq1[0], x);
        //seq2 extended to left of seq1
        } else if( y > 0 && x == 0){
            pre_length = y;
            memcpy(&_merged_seq[0], &_seq2[0], y);
        } else if (x == 0 && y == 0){
            //seq1 & seq2 match on left boundary, do nothging
            pre_length = 0;
        } else {
            printf("Error: Unmatched leading sequence on both contigs.\n");
            exit(1);
        }

        //2. assemble merge during match
        //give seq1 priority (arbitraily) for SNPs, fill in InDels.
        //improvement - use whichever base has higher score for SNPs.
        for(int i=0; i<new_length; ++i){
            if(_gapped_seq1[i] != '-'){
                _merged_seq[i+pre_length] = _gapped_seq1[i];
            } else if (_gapped_seq2[i] != '-'){
                _merged_seq[i+pre_length] = _gapped_seq2[i];
            } else {
                printf("Error: Overlaping Deletiion detected in both aligned sequences.\n");
                exit(1);
            }
        }

        //3. Assemble merge after match
        int post_length;
        //seq1 extends to right of seq2
        if( width > (max_w + 1) ){
            post_length = width-max_w;
            memcpy(&_merged_seq[pre_length + new_length], &_seq1[max_w], post_length);

        //seq2 extends to right of seq1
        } else if( height > (max_h + 1) ){
            post_length = height-max_h;
            memcpy(&_merged_seq[pre_length + new_length], &_seq2[max_h], post_length);

        //seq1 and seq2 both aligned on right boundary
        } else if ( width == max_w && height == max_h ){
            post_length = 0;
        } else {
            printf("Error, unmatched trailing sequence on both contigs.\n");
            exit(1);
        }
        _merged_seq[pre_length + new_length + post_length] = '\0';

        if(DEBUGGING2){
            printf("Merged Sequence:\n%s\n", _merged_seq);
        }
    }

    ~SWMatrix(){
        free(_seq1);
        free(_seq2);
        free(_gapped_seq1);
        free(_gapped_seq2);
        free(_merged_seq);

        for(int i=0; i<height; ++i){
            delete[] h[i];
        }
        delete[] h;
    }

};

#endif
