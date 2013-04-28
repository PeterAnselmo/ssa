
//display debugging info of various verbosity
const bool DEBUGGING = true; 
const bool DEBUGGING2 = false;
const bool DEBUGGING3 = false;

//number of bases in common with edges of consensus & read
//during perfect read assembly stage
const int MIN_OVERLAP = 25;

//when trimming contigs, bases with lower than or eqaual to this quality
//will be removed from edges. Contigs consisting of only bases below
//this quality will be ommitted from contig assembly
const unsigned int CONTIG_TRIM_QUALITY = 2;

//mininum Conitig SW score to consider a match,
//this will need to be adjusted after adjusting match/mismatch scores
const int CONTIG_MATCH_THRESHOLD = 20;

//SW Comparison Scores
const int SW_W_MATCH = 2;
const int SW_W_MISMATCH = -4;
