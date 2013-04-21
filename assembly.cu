#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <algorithm>
#include <iomanip>
#include <list>
#include "contig.cu"
#include "fastqfile.cu"
#include "read.cu"
#include "sw_matrix.cpp"

using namespace std;

const bool DEBUGGING = true;

//number of bases in common with edges of consensus & read
//during perfect read assembly stage
const int MIN_OVERLAP = 20;

//max number of initial perfect match contigs to assemble
const unsigned int CONTIG_CAP = 1000;

//when trimming contigs, bases with lower than or eqaual to this quality
//will be removed from edges. Contigs consisting of only bases below
//this quality will be ommitted from contig assembly
const unsigned int CONTIG_TRIM_QUALITY = 1;

//mininum Conitig SW score to consider a match,
//this will need to be adjusted after adjusting match/mismatch scores
const int CONTIG_MATCH_THRESHOLD = 30;

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
                if( !read->assembled()  && read->seq().find('N') == string::npos){
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

                    cout << "Considering read: " << read->seq() << endl;

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

    void assemble_perfect_contigs_cuda(){
        unsigned int pass = 0;
        //We'll create CONTIG_CAP constigs consisting of perfect matches between reads
        while(pass < CONTIG_CAP){
            Contig c(pass++);

            //find the first unmapped read with no "n" bases to seed the contig
            bool all_mapped = true;
            vector<Read>::iterator read;
            for(read = reads.begin(); read != reads.end(); ++read){
                if( !read->assembled()  && read->seq().find('N') == string::npos){
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

                    cout << "Considering read: " << read->seq() << endl;

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
                if(DEBUGGING3){
                    m.print_matrix();
                }
                if(m.score() > CONTIG_MATCH_THRESHOLD){
                    m.gap_seqs();
                    c1->merge(*c2);
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
        list<Contig>::iterator iter = contigs.begin();
        while(iter != contigs.end()){
            iter->trim(CONTIG_TRIM_QUALITY, reads);
            if(iter->size() == 0){
                if(DEBUGGING){
                    cout << "Deleting low-quality contig: " << iter->id() << endl;
                }
                vector<Read>::iterator read;
                for(read = reads.begin(); read != reads.end(); ++read){
                    if( read->assembled() && read->assem_contig == iter->id() ){
                        read->unassemble();
                    }
                }
                contigs.erase(iter++);
            } else {
                ++iter;
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
        cout << "Assembled " << num_assembled << "/" << reads.size() << " reads ("
             << static_cast<double>(num_assembled) * 100 /reads.size() << "%)" << endl;

    }

    void print_contigs(bool show_reads = true){

        //sort the reads before displaying below contigs
        //reads.sort([](const Read &r1, const Read &r2){ return r1.assem_pos < r2.assem_pos; });

        for(list<Contig>::iterator c = contigs.begin(); c != contigs.end(); ++c){
            cout << "Assembled Contig " << c->id() << ":\n" << c->seq() << endl;
            cout << c->qual() << endl;
            if(show_reads){
                for(vector<Read>::iterator read = reads.begin(); read != reads.end(); ++read){
                    if(read->assembled() && read->assem_contig == c->id()){
                        for(int i=0; i<read->assem_pos; ++i){
                            cout << " ";
                        }
                        cout << read->gapped_seq << endl;
                    }
                }
            }
        }
    }

private:
    void assemble_read(Contig &c, Read &read, unsigned int pos){

        unsigned int overlap_size = min(read.size(), c.size()-pos);

        if(DEBUGGING){
            cout << "Assembling Read: " << read.seq() << " to contig " << c.id() << " at " << pos << endl;
        }
        read.assemble(c.id(), pos);

        SWMatrix sw(c.substr(pos, read.size()), read.seq());
        sw.gap_seqs();
        read.gapped_seq = sw.get_gapped_seq2();

        //if the read extends to the right of the c.seq
        if( overlap_size < read.gapped_seq.size() ){
            if(DEBUGGING){
                cout << "overlap size: " << overlap_size << endl;
                cout << "Read size: " << read.seq().size() << endl;
            }
            for(unsigned int i=0; i<overlap_size; ++i){
                if(c.seq()[pos+i] == read.gapped_seq[i] ){
                    c.inc_qual(pos+i);
                }
            }
            string new_seq = read.gapped_seq.substr(overlap_size);
            if(DEBUGGING){
                cout << "adding " << new_seq << " to reference.\n";
            }
            c.append(new_seq);
        }
    }

    void assemble_perfect_read(Contig &c, Read &read, unsigned int pos){
        unsigned int overlap_size = min(read.size(), c.size()-pos);

        if(DEBUGGING){
            cout << "Assembling Read: " << read.seq() << " to contig " << c.id() << " at " << pos << " "
                 << "overlap size: " << overlap_size << endl;
        }

        //assemble the read here
        read.assemble(c.id(), pos);

        //increment the quality for all of the overlapping bases at the right of the contig
        for(unsigned int j=0; j<overlap_size; ++j){
            c.inc_qual(pos+j);
        }

        if( overlap_size < read.size() ){
            string new_seq = read.seq().substr(overlap_size);

            if(DEBUGGING){
                cout << "adding " << new_seq << " to reference, overlap size: " << overlap_size << ".\n";
            }
            c.append(new_seq);

            if(DEBUGGING){
                cout << "New Reference:\n" << c.seq() << endl;
            }
        }
    }

    void assemble_perfect_read_left(Contig &c, Read &read, unsigned int read_pos){
        unsigned int overlap_size = min(c.size(), read.size()-read_pos);

        if(DEBUGGING){
            cout << "Assembling Read: " << read.seq() << " to contig " << c.id() << " at read pos " << read_pos << " "
                 << "overlap size: " << overlap_size << endl;
        }

        //increment the quality for all of the overlapping bases at the left of the contig
        for(unsigned int j=0; j<overlap_size; ++j){
            c.inc_qual(j);
        }

        if( overlap_size < read.size() ){
            string new_seq = read.seq().substr(0,read_pos);

            if(DEBUGGING){
                cout << "adding " << new_seq << " to beginning reference, overlap size: " << overlap_size << ".\n";
            }
            c.prepend(new_seq);
            c.unshift_aligned_reads(new_seq.size(), reads);

            if(DEBUGGING){
                cout << "New Reference:\n" << c.seq() << endl;
            }
        }
       
        //assemble the read here
        //do it after prepending to sequence, so it is not shifted
        read.assemble(c.id(), 0);
    }

};

#endif
