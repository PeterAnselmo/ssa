#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <algorithm>
#include <iomanip>
#include <list>
#include "contig.cpp"
#include "fastqfile.cpp"
#include "read.cpp"
#include "sw_matrix.cpp"

using namespace std;

//max number of initial perfect match contigs to assemble
const unsigned int CONTIG_CAP = 500;

class Assembly {
public:
    vector<Read> reads;
    list<Contig> contigs;

public:
    Assembly(FastqFile &input_file){
        reads = input_file.reads();
        if(reads.empty()){
            cout << "[ERROR] no reads to align." << endl;
            exit(1);
        }
    }

    //phase 1 - assemble contigs based on perfect overlap in reads - no mismatch allowed
    void assemble_perfect_contigs(){
        unsigned int pass = 0;
        //We'll create CONTIG_CAP constigs consisting of perfect matches between reads
        while(pass < CONTIG_CAP){
            Contig c(pass++);

            //find the first unmapped read with no "n" bases to seed the contig
            bool all_mapped = true;

            vector<Read>::iterator read;
            for(read = reads.begin(); read != reads.end(); ++read){
                if( !read->assembled() && read->find('N') == read->size()){
                    read->assemble(c.id(), 0);
                    c.set_seq(read->seq());
                    all_mapped = false;
                    break;
                }
            }

            //if all the reads were mapped before hitting the contig cap, exit
            if( all_mapped ){
                return;
            }

            if(DEBUGGING) {
                cout << "Starting new contig with sequence:\n" << c.seq() << endl;
            }

            //if any reads were mapped in the last iteration
            bool mapped_read = true;
            while(mapped_read){
                mapped_read = false;

                if(DEBUGGING) {
                    cout << "Restarting at the beginning of read list\n";
                }

                //consider all reads for contig assembly, take the first that matches with > MIN_OVERLAP
                vector<Read>::iterator read;
                for(read = reads.begin(); read != reads.end(); ++read){

                    //if this read was already mapped
                    if( read->assembled() ){
                        continue;
                    }

                    if(DEBUGGING2){
                        printf("Considering read: %s\n", read->seq());
                    }

                    //unsigned int start_pos = c.size() - read->size();
                    unsigned int end_pos = c.size() - MIN_OVERLAP;

                    //compare right side of contig to left side of read
                    //start position depends on whether trying to align reads
                    //for(unsigned int i=start_pos; i<end_pos; ++i){
                    for(unsigned int i=0; i<end_pos; ++i){
                        
                        unsigned int compare_size = min(read->size(), c.size()-i);

                        char *read_substr = read->substr(0,compare_size);
                        char *contig_substr = c.substr(i,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", contig_substr, read_substr);
                        }
                        if( !read->assembled() && strcmp(contig_substr, read_substr) == 0){
                            assemble_perfect_read(c, *read, i);
                            mapped_read = true;
                        }
                        free(read_substr);
                        char *read_rev_substr = read->rev_substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", contig_substr, read_rev_substr);
                        }
                        if( !read->assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                            read->set_rev_comp();
                            assemble_perfect_read(c, *read, i);
                            mapped_read = true;
                        }
                        free(read_rev_substr);
                        free(contig_substr);
                    }
                    
                    //compare left side of contig to right side of red
                    end_pos = read->size() - MIN_OVERLAP;
                    for(unsigned int i=0; i<end_pos; ++i){
                        
                        unsigned int compare_size = min(c.size(), read->size()-i);

                        char *read_substr = read->substr(i,compare_size);
                        char *contig_substr = c.substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", read_substr, contig_substr);
                        }
                        if( !read->assembled() && strcmp(read_substr, contig_substr) == 0 ){
                            assemble_perfect_read_left(c, *read, i);
                            mapped_read = true;
                        }
                        free(read_substr);
                        char *read_rev_substr = read->rev_substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", read_rev_substr, contig_substr);
                        }
                        if( !read->assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                            read->set_rev_comp();
                            assemble_perfect_read_left(c, *read, i);
                            mapped_read = true;
                        }
                        free(read_rev_substr);
                        free(contig_substr);
                    }
                }
            }
            contigs.push_back(c);
        }
    }

    /*
    void assemble_perfect_contigs_cuda(){
        //convert reads to array
        Read *h_reads;
        Read *d_reads;
        h_reads = &reads[0];

        int threads_per_block = 256;
        int num_blocks = reads.size() / threads_per_block + 1;

        //copy congig to device
        //copy reads to device

        unsigned int pass = 0;
        //We'll create CONTIG_CAP constigs consisting of perfect matches between reads
        while(pass < CONTIG_CAP){
            Contig c(pass++);

            //find the first unmapped read with no "n" bases to seed the contig
            bool all_mapped = true;
            vector<Read>::iterator read;
            for(read = reads.begin(); read != reads.end(); ++read){
                if( !read->assembled()  && read->find('N') == read->size()){
                    read->assemble(c.id(), 0);
                    c.seq(read->seq());
                    all_mapped = false;
                    break;
                }
            }

            //if all the reads were mapped before hitting the contig cap, exit
            if( all_mapped ){
                return;
            }


            if(DEBUGGING) {
                cout << "Starting new contig with sequence:\n" << c.set_seq() << endl;
            }

            //if any reads were mapped in the last iteration
            bool mapped_read = true;
            while(mapped_read){
                mapped_read = false;

                if(DEBUGGING) {
                    cout << "Restarting at the beginning of read list\n";
                }

                //compute_overlaps<<<blocks_per_device, threads_per_block>>>(d_reads);

                //consider all reads for contig assembly, take the first that matches with > MIN_OVERLAP
                vector<Read>::iterator read;
                for(read = reads.begin(); read != reads.end(); ++read){

                    //compare left side of contig to right side of red
                    end_pos = read->size() - MIN_OVERLAP;
                    for(unsigned int i=0; i<end_pos; ++i){
                        
                        unsigned int compare_size = min(c.size(), read->size()-i);

                        if(DEBUGGING2){
                            cout << "Considering overlap: " << read->substr(i,compare_size) << "|" << c.substr(0,compare_size) << endl;
                        }
                        if( read->substr(i,compare_size) == c.substr(0,compare_size) ){
                            assemble_perfect_read_left(c, *read, i);
                            mapped_read = true;
                            break;
                        }
                        if( c.substr(i,compare_size) == read->rev_comp().substr(0,compare_size) ){
                            read->set_rev_comp();
                            assemble_perfect_read_left(c, *read, i);
                            mapped_read = true;
                            break;
                        }
                    }
                }
            }
            contigs.push_back(c);
        }
    }
    */

    //phase 2 - assemble contigs to eeach other, allowing mismatches
    void assemble_contigs(){

        //loop over all the contigs, compare them to each other
        list<Contig>::iterator c1;
        for(c1 = contigs.begin(); c1 != contigs.end(); ++c1){

            list<Contig>::iterator c2;
            for(c2 = contigs.begin(); c2 != contigs.end(); ++c2){
                if(c1->id() == c2->id()){
                    continue;
                }
                SWMatrix m(c1->seq(), c2->seq());
                if(DEBUGGING2){
                    printf("Matrix Score: %d\n", m.score());
                }
                if(m.score() >= CONTIG_MATCH_THRESHOLD){
                    m.gap_seqs();
                    if(DEBUGGING){
                        printf("Merged Contigs %d and %d:\n%s\n", c1->id(), c2->id(), m.merged_seq());
                    }
                }
            }


            /*
            //if any reads were mapped in the last iteration
            bool mapped_read = true;
            while(mapped_read){
            mapped_read = false;

                if(DEBUGGING) {
                    cout << "Restarting at the beginning of read list\n";
                }

                //loop over all unmatched reads
                vector<Read>::iterator iter;
                for( iter = reads.begin(); iter != reads.end(); ++iter){
                    Read read = *iter;

                    //if this read was already mapped
                    if( read.assembled() ){
                        continue;
                    }

                    int high_score = 0;
                    int high_pos = 0;
                    vector<Read>::iterator high_iter;

                    unsigned int end_pos = c.size() - MIN_OVERLAP;
                    for(unsigned int i=0; i<end_pos; ++i){
                        unsigned int compare_size = min({read.size(), c.size()-i});

                        SWMatrix sw(c.substr(i,compare_size), read.substr(0,compare_size));
                        SWMatrix rev_sw(c.substr(i,compare_size), read.rev_comp().substr(0,compare_size));

                        if(DEBUGGING2){
                            cout << "Score at pos: " << i << ":" << sw.score() << endl;
                            cout << "RC Score at pos: " << i << ":" << rev_sw.score() << endl;
                        }
                        if( sw.score() > MATCH_THRESHOLD && sw.score() > high_score ){
                            high_score = sw.score();
                            high_iter = iter;
                            high_pos = i;
                        }
                        if( rev_sw.score() > MATCH_THRESHOLD && rev_sw.score() > high_score ){
                            high_score = rev_sw.score();
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
            */
        }
    }

    void trim_contigs(){
        list<Contig>::iterator contig = contigs.begin();
        while(contig != contigs.end()){
            contig->trim(CONTIG_TRIM_QUALITY, reads);
            if(contig->size() == 0){
                if(DEBUGGING){
                    printf("Deleting low-quality contig: %d\n", contig->id());
                }
                contigs.erase(contig++);
            } else {
                ++contig;
            }
        }
    }

    void print_report(){
        unsigned int num_assembled = 0;
        for(vector<Read>::iterator read = reads.begin(); read != reads.end(); ++read){
            if( read->assembled() ){
                ++num_assembled;
            }
        }
        printf("Assembled %u/%u reads(%f%%)", num_assembled, 
                                              static_cast<unsigned int>(reads.size()), 
                                              static_cast<double>(num_assembled) * 100 /reads.size());

    }

    void print_contigs(bool show_reads = true){

        if(show_reads){
            //sort the reads before displaying below contigs
            sort(reads.begin(), reads.end());
        }

        for(list<Contig>::iterator c = contigs.begin(); c != contigs.end(); ++c){
            printf("Assembled Contig %d:\n%s\n%s\n", c->id(), c->seq(), c->qual());
            if(show_reads){
                for(vector<Read>::iterator read = reads.begin(); read != reads.end(); ++read){
                    if(read->assembled() && read->contig() == c->id()){
                        for(unsigned int i=0; i<read->position(); ++i){
                            printf(" ");
                        }
                        printf("%s\n",read->gapped_seq());
                    }
                }
            }
        }
    }

private:
    /*
    void assemble_read(Contig &c, Read &read, unsigned int pos){

        unsigned int overlap_size = min(read.size(), c.size()-pos);

        if(DEBUGGING){
            cout << "Assembling Read: " << read.seq() << " to contig " << c.id() << " at " << pos << endl;
        }
        read.assemble(c.id(), pos);

        SWMatrix sw(c.substr(pos, read.size()), read.seq());
        sw.gap_seqs();
        strcpy(read.gapped_seq, sw.get_gapped_seq2().c_str());

        //if the read extends to the right of the c.seq
        if( overlap_size < read.gapped_size() ){
            if(DEBUGGING){
                cout << "overlap size: " << overlap_size << endl;
                cout << "Read size: " << read.size() << endl;
            }
            for(unsigned int i=0; i<overlap_size; ++i){
                if(c.seq()[pos+i] == read.gapped_seq[i] ){
                    c.inc_qual(pos+i);
                }
            }
            char* new_seq = read.gapped_substr(overlap_size);
            if(DEBUGGING){
                cout << "adding " << new_seq << " to reference.\n";
            }
            c.append(new_seq);
        }
    }
    */

    void assemble_perfect_read(Contig &c, Read &read, unsigned int pos){
        unsigned int overlap_size = min(read.size(), c.size()-pos);

        if(DEBUGGING){
            printf("Assembling read: %s to contig %d at %d, overlap size: %u\n", read.seq(), c.id(), pos, overlap_size);
        }

        //assemble the read here
        read.assemble(c.id(), pos);

        //increment the quality for all of the overlapping bases at the right of the contig
        for(unsigned int j=0; j<overlap_size; ++j){
            c.inc_qual(pos+j);
        }

        if( overlap_size < read.size() ){
            char *new_seq = read.substr(overlap_size);

            if(DEBUGGING){
                printf("Adding %s to end of reference\n", new_seq);
            }
            c.append(new_seq);
            free(new_seq);

            if(DEBUGGING){
                printf("New Reference:\n%s\n", c.seq());
            }
        }
    }

    void assemble_perfect_read_left(Contig &c, Read &read, unsigned int read_pos){
        unsigned int overlap_size = min(c.size(), read.size()-read_pos);

        if(DEBUGGING){
            printf("Assembling Read: %s to contig %d at pos %d with overlap %d\n", read.seq(), c.id(), read_pos, overlap_size);
        }

        //increment the quality for all of the overlapping bases at the left of the contig
        for(unsigned int j=0; j<overlap_size; ++j){
            c.inc_qual(j);
        }

        if( overlap_size < read.size() ){
            char *new_seq = read.substr(0,read_pos);

            if(DEBUGGING){
                printf("Adding %s to beginning of reference\n", new_seq);
            }
            c.prepend(new_seq);
            c.unshift_aligned_reads(strlen(new_seq), reads);

            if(DEBUGGING){
                printf("New Reference:\n%s\n", c.seq());
            }
            free(new_seq);
        }
       
        //assemble the read here
        //do it after prepending to sequence, so it is not shifted
        read.assemble(c.id(), 0);
    }

};

/*
__global__ void compute_overlaps(Read *reads, int num_reads, int contig_id, int contig_size ){
    int tid = blockIdx.x*blockDim.x + threadIdx.x;

    //there may be a few extra threads called, make sure in range
    if( tid < num_reads ){
        Read read = reads[tid];

        //if this read was already mapped
        if( read.assembled() ){
            return;
        }

        unsigned int end_pos = c.size() - MIN_OVERLAP;

        //compare right side of contig to left side of read
        for(unsigned int i=0; i<end_pos; ++i){
            
            unsigned int compare_size = min(read->size(), c.size()-i);

            if(DEBUGGING2){
                cout << "Considering overlap: " << c.substr(i,compare_size) << "|" << read->substr(0,compare_size) << endl;
            }
            if( c.substr(i,compare_size) == read->substr(0,compare_size) ){
                assemble_perfect_read(c, *read, i);
                mapped_read = true;
                break;
            }
            if( c.substr(i,compare_size) == read->rev_comp().substr(0,compare_size) ){
                read->set_rev_comp();
                assemble_perfect_read(c, *read, i);
                mapped_read = true;
                break;
            }
        }


    }
}
*/


#endif
