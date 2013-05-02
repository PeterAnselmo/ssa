#ifndef SW_MATRIX_H
#define SW_MATRIX_H

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "settings.cpp"

using namespace std;


class SWMatrix {
private:
    //matrix to contain scores.
    int **h;

    Contig* c1;
    Contig* c2;

    //size of the matrix.  Width will be the size of seq1 & height will be size of seq2
    int height;
    int width;

    //coordinate of the highest score in matrix (only last row/column considered (Global matching))
    int max_h, max_w; 

    //These are the overlapping portions of the sequence that overlap with gaps (-) added after alignment
    char* _gapped_seq1;
    char* _gapped_seq2;

    //this is the flattened version of the gapped sequences, with all bases filled, and 
    //SNPs chosen based on which sequence had the higher quality
    char* _merged_seq;
    char* _merged_qual;

    //this is the merged sequence with flanking sequence before and after the match
    //from the two reads combined into one long read.
    char* _complete_seq;
    char* _complete_qual;


public:
    SWMatrix(Contig &cont1, Contig &cont2){
        c1 = &cont1;
        c2 = &cont2;

        //coincidentally, the matrix will have the same dimmensions as the allocated 
        //strings, since there is an extra zero row and a null character respectively.
        width = c1->size() + 1;
        height = c2->size() + 1;
        
        _gapped_seq1 = NULL;
        _gapped_seq2 = NULL;
        _merged_seq = NULL;
        _merged_qual = NULL;
        _complete_seq = NULL;
        _complete_qual = NULL;
        
        max_h = max_w = 0;

        if(DEBUGGING2){
            printf("Constructing matrix of c%dxc%d (%dx%d)\n", c1->id(), c2->id(), width, height);
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

    char* complete_seq(){
        return _complete_seq;
    }
    char* complete_qual(){
        return _complete_qual;
    }

    void compute_matrix(){
        max_h = max_w = 0;

        int max_val = 0;
        for(int i=1; i<height; i++){
            for(int j=1; j<width; ++j){
                int w = (c1->seq()[j-1] == c2->seq()[i-1]) ? SW_W_MATCH : SW_W_MISMATCH;
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
            cout << setw(padding+1) << c1->seq()[i];
        }
        cout << endl;
        cout << setw(padding-1) << '-';
        for(int i=0; i<width; ++i){
            cout << setw(padding+1) << '0';
        }
        cout << endl;
        for(int i=1; i<height; ++i){
            cout << setw(padding) << c2->seq()[i-1];
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
    void merge_seqs(){
        //we'll allocate the worst case scenario (length of both strings) to avoid the need to reallocate
        _gapped_seq1 = (char*)malloc(max_w+max_h+1);
        _gapped_seq2 = (char*)malloc(max_w+max_h+1);
        _merged_seq = (char*)malloc(max_w+max_h+1);
        _merged_qual = (char*)malloc(max_w+max_h+1);
        _complete_seq = (char*)malloc(width+height+1);
        _complete_qual = (char*)malloc(width+height+1);

        //this will be the position of the gapped sequence as we assign
        int current_position = max_w+max_h;
        
        //x & y will trace the path from lower right to top left
        //currently gives preferences as follows: match/mismatch, deletion, insertion
        //an improvment would be to use the one with the highest prev score
        int x = max_w;
        int y = max_h;

        if( DEBUGGING2){
            printf("Computing Gapped Sequences, starting with %d,%d...\n", y,x);
        }

        while(x > 0 && y > 0){

            //diagonal - match
            //note, string position is -1 from grid position
            if( (c1->seq()[x-1] == c2->seq()[y-1]) && (h[y][x] == (h[y-1][x-1] + SW_W_MATCH)) ){
                --x;
                --y;

                _gapped_seq1[current_position] = c1->seq()[x];
                _gapped_seq2[current_position] = c2->seq()[y];
                _merged_seq[current_position] = c1->seq()[x];

                int temp_qual = c1->qual()[x] + c2->qual()[y] - QUAL_OFFSET;
                if(temp_qual > MAX_QUAL){ temp_qual = MAX_QUAL; }
                _merged_qual[current_position] = temp_qual;

                if(DEBUGGING3){ printf("Traceback: Diagonal Match, %cvs%c, %cvs%c, final qual: %c\n",
                        c1->seq()[x], 
                        c2->seq()[y], 
                        c1->qual()[x], 
                        c2->qual()[y], 
                        _merged_qual[current_position]); }


            //diagonal - mismatch       
            } else if( (c1->seq()[x-1] != c2->seq()[y-1]) && (h[y][x] == (h[y-1][x-1] + SW_W_MISMATCH))) {
                --x;
                --y;

                _gapped_seq1[current_position] = c1->seq()[x];
                _gapped_seq2[current_position] = c2->seq()[y];
                if( c1->qual()[x] >= c2->qual()[x] ){
                    _merged_seq[current_position] = c1->seq()[x];
                    _merged_qual[current_position] = c1->qual()[x];
                } else {
                    _merged_seq[current_position] = c2->seq()[y];
                    _merged_qual[current_position] = c2->qual()[y];
                }

                if(DEBUGGING3){ printf("Traceback: Diagonal Mismatch, %cvs%c, %cvs%c, final qual: %c\n",
                        c1->seq()[x], 
                        c2->seq()[y], 
                        c1->qual()[x], 
                        c2->qual()[y], 
                        _merged_qual[current_position]); }

            //moving left - deletion in second sequence
            } else if ( h[y][x] == (h[y][x-1] + SW_W_MISMATCH) ){
                --x;
                
                if(DEBUGGING3){ printf("Traceback: Left\n");}

                _gapped_seq1[current_position] = c1->seq()[x];
                _gapped_seq2[current_position] = '-';
                _merged_seq[current_position] = c1->seq()[x];
                _merged_qual[current_position] = c1->qual()[x];

            //moving up - insertion in second sequence
            } else if ( h[y][x] == (h[y-1][x] + SW_W_MISMATCH) ) {
                --y;

                if(DEBUGGING3){ printf("Traceback: Up\n"); }

                _gapped_seq1[current_position] = '-';
                _gapped_seq2[current_position] = c2->seq()[y];
                _merged_seq[current_position] = c2->seq()[y];
                _merged_qual[current_position] = c2->qual()[y];

            } else if( h[y][x] == 0 ){
                --x;
                --y;
                
                if(DEBUGGING3){ printf("Traceback Lost Path: Guessing Diagonal\n"); }

                _gapped_seq1[current_position] = c1->seq()[x];
                _gapped_seq2[current_position] = c2->seq()[y];
                _merged_seq[current_position] = c1->seq()[x];
                _merged_qual[current_position] = '!';
                
            } else {
                printf("Error, could not compute valid traceback");
                exit(1);
            }
            --current_position;
        }
        //Fencepost
        _gapped_seq1[current_position] = c1->seq()[x];
        _gapped_seq2[current_position] = c2->seq()[y];
        if(c1->seq()[x] == c2->seq()[y]){
            _merged_seq[current_position] = c1->seq()[x];
            _merged_qual[current_position] = c1->qual()[x] + c2->qual()[y] - QUAL_OFFSET;
        } else {
            if( c1->qual()[x] >= c2->qual()[x] ){
                _merged_seq[current_position] = c1->seq()[x];
                _merged_qual[current_position] = c1->qual()[x];
            } else {
                _merged_seq[current_position] = c2->seq()[y];
                _merged_qual[current_position] = c2->qual()[y];
            }
        }

        //move all backtracked sequences to start of their array
        int new_length = max_w + max_h - current_position;
        ++current_position;
        memmove(&_gapped_seq1[0], &_gapped_seq1[current_position], new_length);
        _gapped_seq1[new_length] = '\0';
        memmove(&_gapped_seq2[0], &_gapped_seq2[current_position], new_length);
        _gapped_seq2[new_length] = '\0';
        memmove(&_merged_seq[0], &_merged_seq[current_position], new_length);
        _merged_seq[new_length] = '\0';
        memmove(&_merged_qual[0], &_merged_qual[current_position], new_length);
        _merged_qual[new_length] = '\0';

        if(DEBUGGING3){
            printf("Gapped Contigs:\n%s\n%s\n", _gapped_seq1, _gapped_seq2);
            printf("Merged Result:\n%s\n%s\n", _merged_seq, _merged_qual);
        }

        //1. assemble merge before match
        int pre_length;
        //seq1 extended to left of seq2
        if(x > 0 && y == 0){
            pre_length = x;
            memcpy(&_complete_seq[0], &c1->seq()[0], x);
            memcpy(&_complete_qual[0], &c1->qual()[0], x);

        //seq2 extended to left of seq1
        } else if( y > 0 && x == 0){
            pre_length = y;
            memcpy(&_complete_seq[0], &c2->seq()[0], y);
            memcpy(&_complete_qual[0], &c2->qual()[0], y);

        //seq1 & seq2 match on left boundary, do nothing
        } else if (x == 0 && y == 0){
            pre_length = 0;
        } else {
            printf("Error: Unmatched leading sequence on both contigs.\n");
            exit(1);
        }

        //2. assemble merge during match
        memcpy(&_complete_seq[pre_length], &_merged_seq[0], new_length);
        memcpy(&_complete_qual[pre_length], &_merged_qual[0], new_length);

        //3. Assemble merge after match
        int post_length;
        //seq1 extends to right of seq2
        if( width > (max_w + 1) ){
            post_length = width-max_w;
            memcpy(&_complete_seq[pre_length + new_length], &c1->seq()[max_w], post_length);
            memcpy(&_complete_qual[pre_length + new_length], &c1->qual()[max_w], post_length);

        //seq2 extends to right of seq1
        } else if( height > (max_h + 1) ){
            post_length = height-max_h;
            memcpy(&_complete_seq[pre_length + new_length], &c2->seq()[max_h], post_length);
            memcpy(&_complete_qual[pre_length + new_length], &c2->qual()[max_h], post_length);

        //seq1 and seq2 both aligned on right boundary
        } else if ( width == (max_w+1) && height == (max_h+1) ){
            post_length = 0;
        } else {
            printf("Error, unmatched trailing sequence on both contigs.\n");
            exit(1);
        }
        _complete_seq[pre_length + new_length + post_length] = '\0';
        _complete_qual[pre_length + new_length + post_length] = '\0';

        if(DEBUGGING){
            printf("Merged Sequence:\n%s\n%s\n", _complete_seq, _complete_qual);
        }
    }

    ~SWMatrix(){
        free(_gapped_seq1);
        free(_gapped_seq2);
        free(_merged_seq);
        free(_merged_qual);
        free(_complete_seq);
        free(_complete_qual);

        for(int i=0; i<height; ++i){
            delete[] h[i];
        }
        delete[] h;
    }

};

#endif
