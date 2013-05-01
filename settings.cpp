#ifndef SETTINGS_CPP
#define SETTINGS_CPP

//display debugging info of various verbosity
const bool DEBUGGING = false; 
const bool DEBUGGING2 = false;
const bool DEBUGGING3 = false;

//number of bases in common with edges of consensus & read
//during the first perfect read assembly stage
const int MIN_OVERLAP = 20;

//when trimming contigs, bases with lower than this setting
//will be removed from edges. Contigs consisting of only bases below
//this quality will be ommitted from contig assembly
const unsigned int CONTIG_TRIM_QUALITY = 5;

//max number of initial perfect match contigs to assemble
const int CONTIG_CAP = 500;

//SW Comparison Scores
const int SW_W_MATCH = 2;
const int SW_W_MISMATCH = -4;

//mininum Conitig SW score to consider a match,
//this will need to be adjusted after adjusting match/mismatch scores
const int CONTIG_MATCH_THRESHOLD = 24;

//number of bases per line in output fasta file.
const int FASTA_LINE_WIDTH = 80;

//quality is represented by a single ascii character.  The range of safely printable cahracters
//ranges from 33 (!) to 126(~).  Don't assign qualites out of this range.
///Don't change these.
const int MAX_QUAL = 126;
const int QUAL_OFFSET = 32;

#endif
