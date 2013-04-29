#ifndef SETTINGS_CPP
#define SETTINGS_CPP

//display debugging info of various verbosity
const bool DEBUGGING = true; 
const bool DEBUGGING2 = true;
const bool DEBUGGING3 = true;

//truncate all reads to this many chars
//would be nice to get ride of this.
const int READ_SIZE = 30;
//conversion from quality of 1 to lowest printable ascii char (33-!)
//Don't change this.
const int QUAL_OFFSET = 32;

//number of bases in common with edges of consensus & read
//during perfect read assembly stage
const int MIN_OVERLAP = 20;

//when trimming contigs, bases with lower than this setting
//will be removed from edges. Contigs consisting of only bases below
//this quality will be ommitted from contig assembly
const unsigned int CONTIG_TRIM_QUALITY = 2;

//max number of initial perfect match contigs to assemble
const unsigned int CONTIG_CAP = 500;

//SW Comparison Scores
const int SW_W_MATCH = 2;
const int SW_W_MISMATCH = -4;

//mininum Conitig SW score to consider a match,
//this will need to be adjusted after adjusting match/mismatch scores
const int CONTIG_MATCH_THRESHOLD = 24;

#endif
