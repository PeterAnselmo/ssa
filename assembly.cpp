#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <algorithm>
#include <iomanip>
#include <list>
#include <omp.h>
#include "contig.cpp"
#include "fastq.cpp"
#include "read.cpp"
#include "sw_matrix.cpp"

using namespace std;

class Assembly {
public:
    vector<Read> reads;
    list<Contig> contigs;

public:
    Assembly(vector<Read> &new_reads){
        reads = new_reads;
        if(reads.empty()){
            printf("[ERROR] no reads to align.\n");
            exit(1);
        }
    }

    //phase 1 - assemble contigs based on perfect overlap in reads - no mismatch allowed
    void assemble_perfect_contigs(){

        //We'll create CONTIG_CAP contigs consisting of perfect matches between reads
        //each contig can (for the most part) assemble in parallel
        #pragma omp parallel for
        for(int pass = 0; pass<CONTIG_CAP; ++pass){
            Contig c(pass);

            //whether there are any valid umapped reads to seed a contig
            bool all_mapped = true;

            //Possibility exists for race condition between checking if a read is
            //mapped, and using it to start a contig.  Need mutex here.
            #pragma omp critical(seed_contig)
            {
                for(unsigned int i=pass; i<reads.size(); ++i){
                    if( !reads[i].assembled() && reads[i].find('N') == reads[i].size()){
                        reads[i].assemble(c.id(), 0);
                        c.set_seq(reads[i].seq());
                        all_mapped = false;
                        break;
                    }
                }
            }//end mutex

            //if all the reads were mapped before hitting the contig cap, exit
            if( all_mapped ){
                continue;
            }

            if(DEBUGGING) {
                printf("Starting new contig with sequence: %s\n", c.seq());
            }

            //if any reads were mapped in the last iteration
            bool mapped_read = true;

            //keep looping and adding to contig until no more reads can be mapped
            while(mapped_read){
                mapped_read = false;

                if(DEBUGGING2) {
                    printf("Restarting at the beginning of read list\n");
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
                            #pragma omp critical(assem_perfect)
                            {
                                if( !read->assembled() && strcmp(contig_substr, read_substr) == 0){
                                    assemble_perfect_read(c, *read, i);
                                    mapped_read = true;
                                }
                            }//end mutex
                        }
                        free(read_substr);
                        char *read_rev_substr = read->rev_substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", contig_substr, read_rev_substr);
                        }
                        if( !read->assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                            #pragma omp critical(assem_perfect)
                            {
                                if( !read->assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                                    read->set_rev_comp();
                                    assemble_perfect_read(c, *read, i);
                                    mapped_read = true;
                                }
                            }//end mutex
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
                            #pragma omp critical(assem_perfect)
                            {
                                if( !read->assembled() && strcmp(read_substr, contig_substr) == 0 ){
                                    assemble_perfect_read_left(c, *read, i);
                                    mapped_read = true;
                                }
                            }//end mutex
                        }
                        free(read_substr);
                        char *read_rev_substr = read->rev_substr(0,compare_size);
                        if(DEBUGGING2){
                            printf("Considering overlap: %s | %s\n", read_rev_substr, contig_substr);
                        }
                        if( !read->assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                            #pragma omp critical(assem_perfect)
                            {
                                if( !read->assembled() && strcmp(contig_substr, read_rev_substr) == 0 ){
                                    read->set_rev_comp();
                                    assemble_perfect_read_left(c, *read, i);
                                    mapped_read = true;
                                }
                            }//end mutex
                        }
                        free(read_rev_substr);
                        free(contig_substr);
                    }
                }
            }
            #pragma omp critical(push_contig)
            contigs.push_back(c);
        }
    }


    //phase 2 - assemble contigs to eeach other, allowing mismatches
    void assemble_contigs(){
        //keep looping as long as there a high quality matches
        bool had_merge = true;

        while(had_merge){
            had_merge = false;
            if(DEBUGGING){
                printf("Restarting Contig Merge Loop.\n");
            }
            //
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

    //The end result of assembly. Ideally, there will be only one contig remaining after merging
    //if the assembly was not awesome, there can be several contigs, in which case, we take the longest
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
            contig->trim(CONTIG_TRIM_QUALITY, reads);
            if(contig->size() == 0){
                if(DEBUGGING){
                    printf("Deleting low-quality contig: %d\n", contig->id());
                }
                //this will erase the contig while at the same time increment the iterator
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
        printf("Assembled %u/%u reads(%f%%)\n", num_assembled, 
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


#endif
