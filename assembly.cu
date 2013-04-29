#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <algorithm>
#include <iomanip>
#include <list>
#include "stdio.h"
#include "contig.cu"
#include "fastqfile.cu"
#include "read.cu"
#include "sw_matrix.cu"

using namespace std;

__host__ __device__ void assemble_perfect_read(Contig *c, Read &read, unsigned int pos);
__host__ __device__ void assemble_perfect_read_left(Contig *c, Read &read, unsigned int pos);
__global__ void map_reads_on_device(Contig *c, Read *reads, int num_reads);

class Assembly {
public:
    Read* reads;
    int num_reads;
    list<Contig> contigs;

public:
    Assembly(FastqFile &input_file){
        reads = input_file.reads();
        num_reads = input_file.num_reads();
        if(num_reads == 0){
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
            Contig *cPtr = &c;

            //find the first unmapped read with no "n" bases to seed the contig
            bool all_mapped = true;

            for(int i=0; i<num_reads; ++i){
                if( !reads[i].assembled() && reads[i].find('N') == reads[i].size()){
                    reads[i].assemble(c.id(), 0);
                    c.set_seq(reads[i].seq());
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
                for(int i=0; i<num_reads; ++i){

                    //if this read was already mapped
                    if( reads[i].assembled() ){
                        continue;
                    }

                    if(DEBUGGING2){
                        printf("Considering read: %s\n", reads[i].seq());
                    }

                    //unsigned int start_pos = c.size() - reads[i].size();
                    unsigned int end_pos = c.size() - MIN_OVERLAP;

                    //compare right side of contig to left side of read
                    //start position depends on whether trying to align reads
                    //for(unsigned int i=start_pos; i<end_pos; ++i){
                    for(unsigned int j=0; j<end_pos; ++j){
                        
                        unsigned int compare_size = min(reads[i].size(), c.size()-j);

                        char *read_substr = reads[i].substr(0,compare_size);
                        char *contig_substr = c.substr(j,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", contig_substr, read_substr);
                        }
                        if( !reads[i].assembled() && strcmp(contig_substr, read_substr) == 0){
                            assemble_perfect_read(cPtr, reads[i], j);
                            mapped_read = true;
                        }
                        free(read_substr);
                        char *read_rev_substr = reads[i].rev_substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", contig_substr, read_rev_substr);
                        }
                        if( !reads[i].assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                            reads[i].set_rev_comp();
                            assemble_perfect_read(cPtr, reads[i], j);
                            mapped_read = true;
                        }
                        free(read_rev_substr);
                        free(contig_substr);
                    }
                    
                    //compare left side of contig to right side of read
                    end_pos = reads[i].size() - MIN_OVERLAP;
                    for(unsigned int j=0; j<end_pos; ++j){
                        
                        unsigned int compare_size = min(c.size(), reads[i].size()-j);

                        char *read_substr = reads[i].substr(j,compare_size);
                        char *contig_substr = c.substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", read_substr, contig_substr);
                        }
                        if( !reads[i].assembled() && strcmp(read_substr, contig_substr) == 0 ){
                            assemble_perfect_read_left(cPtr, reads[i], j);
                            mapped_read = true;
                        }
                        free(read_substr);
                        char *read_rev_substr = reads[i].rev_substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", read_rev_substr, contig_substr);
                        }
                        if( !reads[i].assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                            reads[i].set_rev_comp();
                            assemble_perfect_read_left(cPtr, reads[i], j);
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

    void assemble_perfect_contigs_cuda(){

        //copy all reads to device
        Read *d_reads;
        int reads_size = num_reads * sizeof(Read);
        cudaMalloc( (void**)&d_reads, reads_size);
        cudaMemcpy( d_reads, reads, reads_size, cudaMemcpyHostToDevice);

        int threads_per_block = 256;
        int num_blocks = num_reads / threads_per_block + 1;

        unsigned int pass = 0;
        //We'll create CONTIG_CAP constigs consisting of perfect matches between reads
        while(pass < CONTIG_CAP){
            Contig c(pass++);
            Contig *cPtr = &c;

            //find the first unmapped read with no "n" bases to seed the contig
            bool all_mapped = true;

            for(int i=0; i<num_reads; ++i){
                if( !reads[i].assembled() && reads[i].find('N') == reads[i].size()){
                    reads[i].assemble(c.id(), 0);
                    c.set_seq(reads[i].seq());
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

            //move contig onto device
            Contig *d_contig;
            cudaMalloc( (void**)&d_contig, sizeof(Contig));
            cudaMemcpy( d_contig, &c, sizeof(Contig), cudaMemcpyHostToDevice);

            //if any reads were mapped in the last iteration
            bool mapped_read = true;
            while(mapped_read){
                mapped_read = false;

                if(DEBUGGING) {
                    cout << "Restarting at the beginning of read list\n";
                }

                //consider all reads for contig assembly, take the first that matches with > MIN_OVERLAP
                map_reads_on_device<<<num_blocks, threads_per_block>>>(d_contig, d_reads, num_reads);

                cudaDeviceSynchronize();
                cudaMemcpy( cPtr, d_contig, sizeof(Contig), cudaMemcpyDeviceToHost);

                contigs.push_back(*cPtr);
            }
        }

        //copy reads back from device
        cudaMemcpy( reads, d_reads, reads_size, cudaMemcpyDeviceToHost);
    }

    //phase 2 - assemble contigs to eeach other, allowing mismatches
    void assemble_contigs(){
        //keep looping as long as there a high quality matches
        bool had_merge = true;

        while(had_merge){
            if(DEBUGGING){
                printf("Restarting Contig Merge Loop.\n");
            }
            had_merge = false;
            //loop over all the contigs, compare them to each other
            list<Contig>::iterator c1;
            for(c1 = contigs.begin(); c1 != contigs.end(); ++c1){
                if(DEBUGGING){
                    printf("Comparing contig %d against others\n", c1->id());
                }
                //needs to be while rather then for to handle deleting elements real time
                list<Contig>::iterator c2 = contigs.begin();;
                while(c2 != contigs.end()){
                    if(c1->id() == c2->id()){
                        ++c2;
                        continue;
                    }

                    //compute the score matrix using the sequences from the two contigs
                    SWMatrix m(*c1, *c2);
                    Contig* rev_c2 = c2->rev_comp();
                    SWMatrix m2(*c1, *rev_c2);
                    if(DEBUGGING2){
                        printf("Matrix Score: %d\n", m.score());
                    }

                    //if these two contigs are a match, merge the second one
                    //into the first one and delete the second.
                    if(m.score() >= CONTIG_MATCH_THRESHOLD){
                        if(DEBUGGING){
                            printf("Merging Contigs %d & %d\n", c1->id(), c2->id());
                        }
                        had_merge = true;
                        m.merge_seqs();
                        c1->set_seq(m.complete_seq(), false);
                        c1->set_qual(m.complete_qual());
                        contigs.erase(c2++);

                    //if the rev_comp of the second matches, merge it into the first
                    } else if (m2.score() > CONTIG_MATCH_THRESHOLD ){
                        if(DEBUGGING){
                            printf("Merging Contigs %d & rev %d\n", c1->id(), c2->id());
                        }
                        had_merge = true;
                        m2.merge_seqs();
                        c1->set_seq(m2.complete_seq(), false);
                        c1->set_qual(m2.complete_qual());
                        contigs.erase(c2++);
                    } else {
                        ++c2;
                    }
                    delete(rev_c2);
                }//inner contigs
            }//outer contigs
        }//while had merge
    }

    char* final_seq(){
        list<Contig>::iterator max;
        unsigned int max_size = 0;

        for(list<Contig>::iterator c = contigs.begin(); c != contigs.end(); ++c){
            if(c->size() > max_size){
                max_size = c->size();
                max = c;
            }
        }
        if(max_size == 0){
            printf("Error, unable to select final contig");
        }
        return max->seq();
    }


    void trim_contigs(){
        list<Contig>::iterator contig = contigs.begin();
        while(contig != contigs.end()){
            contig->trim(CONTIG_TRIM_QUALITY, reads, num_reads);
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
        for(int i=0; i<num_reads; ++i){
            if( reads[i].assembled() ){
                ++num_assembled;
            }
        }
        printf("Assembled %u/%u reads(%f%%)\n", num_assembled, 
                                              num_reads,
                                              static_cast<double>(num_assembled) * 100 /num_reads);

    }

    void print_contigs(bool show_reads = true){

        if(show_reads){
            //sort the reads before displaying below contigs
            //sort(reads.begin(), reads.end());
        }

        for(list<Contig>::iterator c = contigs.begin(); c != contigs.end(); ++c){
            printf("Assembled Contig %d:\n%s\n%s\n", c->id(), c->seq(), c->qual());
            if(show_reads){
                for(int i=0; i<num_reads; ++i){
                    if(reads[i].assembled() && reads[i].contig() == c->id()){
                        for(unsigned int j=0; j<reads[i].position(); ++j){
                            printf(" ");
                        }
                        printf("%s\n",reads[i].gapped_seq());
                    }
                }
            }
        }
    }

private:

};

    __host__ __device__ void assemble_perfect_read(Contig *c, Read &read, unsigned int pos){
        unsigned int overlap_size;
        if(read.size() < (c->size()-pos)){
           overlap_size = read.size();
        } else {
           overlap_size = c->size()-pos;
        }

        if(DEBUGGING){
            printf("Assembling read: %s to contig %d at %d, overlap size: %u\n", read.seq(), c->id(), pos, overlap_size);
        }

        //assemble the read here
        read.assemble(c->id(), pos);

        //increment the quality for all of the overlapping bases at the right of the contig
        for(unsigned int j=0; j<overlap_size; ++j){
            c->inc_qual(pos+j);
        }

        if( overlap_size < read.size() ){
            char *new_seq = read.substr(overlap_size);

            if(DEBUGGING){
                printf("Adding %s to end of reference\n", new_seq);
            }
            c->append(new_seq);
            free(new_seq);

            if(DEBUGGING){
                printf("New Reference:\n%s\n", c->seq());
            }
        }
    }

    __host__ __device__ void assemble_perfect_read_left(Contig *c, Read &read, unsigned int pos){
        unsigned int overlap_size;
        if(read.size() < (c->size()-pos)){
           overlap_size = read.size();
        } else {
           overlap_size = c->size()-pos;
        }


        if(DEBUGGING){
            printf("Assembling Read: %s to contig %d at pos %d with overlap %d\n", read.seq(), c->id(), pos, overlap_size);
        }

        //increment the quality for all of the overlapping bases at the left of the contig
        for(unsigned int j=0; j<overlap_size; ++j){
            c->inc_qual(j);
        }

        if( overlap_size < read.size() ){
            char *new_seq = read.substr(0,pos);

            if(DEBUGGING){
                printf("Adding %s to beginning of reference\n", new_seq);
            }
            c->prepend(new_seq);
            #ifndef __CUDA_ARCH__
//                c->unshift_aligned_reads(strlen(new_seq), reads, num_reads);
            #else
 //               c->unshift_aligned_reads(cudaStrlen(new_seq), reads, num_reads);
            #endif

            if(DEBUGGING){
                printf("New Reference:\n%s\n", c->seq());
            }
            free(new_seq);
        }
       
        //assemble the read here
        //do it after prepending to sequence, so it is not shifted
        read.assemble(c->id(), 0);
    }


    __global__ void map_reads_on_device(Contig *c, Read *reads, int num_reads){
        int tid = blockIdx.x*blockDim.x + threadIdx.x;

        //there may be a few extra threads called, make sure in range
        if( tid >= num_reads ){
            return;
        }

        //if this read was already mapped
        if( reads[tid].assembled() ){
            return;
        }

        if(DEBUGGING2){
            printf("Considering read: %s\n", reads[tid].seq());
        }

        //unsigned int start_pos = c.size() - reads[tid].size();
        unsigned int end_pos = c->size() - MIN_OVERLAP;

        //compare right side of contig to left side of read
        //start position depends on whether trying to align reads
        //for(unsigned int i=start_pos; i<end_pos; ++i){
        for(unsigned int j=0; j<end_pos; ++j){
            
            unsigned int compare_size = min(reads[tid].size(), c->size()-j);

            char *read_substr = reads[tid].substr(0,compare_size);
            char *contig_substr = c->substr(j,compare_size);
            if(DEBUGGING2){
                printf("Considering overlap: %s | %s\n", contig_substr, read_substr);
            }
            if( !reads[tid].assembled() && cudaStrcmp(contig_substr, read_substr) == 0){
                assemble_perfect_read(c, reads[tid], j);
            }
            free(read_substr);
            char *read_rev_substr = reads[tid].rev_substr(0,compare_size);
            if(DEBUGGING2){
                printf("Considering overlap: %s | %s\n", contig_substr, read_rev_substr);
            }
            if( !reads[tid].assembled() && cudaStrcmp(contig_substr, read_rev_substr) == 0 ){
                reads[tid].set_rev_comp();
                assemble_perfect_read(c, reads[tid], j);
            }
            free(read_rev_substr);
            free(contig_substr);
        }
        
        //compare left side of contig to right side of read
        end_pos = reads[tid].size() - MIN_OVERLAP;
        for(unsigned int j=0; j<end_pos; ++j){
            
            unsigned int compare_size = min(c->size(), reads[tid].size()-j);

            char *read_substr = reads[tid].substr(j,compare_size);
            char *contig_substr = c->substr(0,compare_size);
            if(DEBUGGING2){
                printf("Considering overlap: %s | %s\n", read_substr, contig_substr);
            }
            if( !reads[tid].assembled() && cudaStrcmp(read_substr, contig_substr) == 0 ){
                assemble_perfect_read_left(c, reads[tid], j);
            }
            free(read_substr);
            char *read_rev_substr = reads[tid].rev_substr(0,compare_size);
            if(DEBUGGING2){
                printf("Considering overlap: %s | %s\n", read_rev_substr, contig_substr);
            }
            if( !reads[tid].assembled() && cudaStrcmp(contig_substr, read_rev_substr) == 0 ){
                reads[tid].set_rev_comp();
                assemble_perfect_read_left(c, reads[tid], j);
            }
            free(read_rev_substr);
            free(contig_substr);
        }
    }


#endif
