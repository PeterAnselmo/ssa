SSA - Super Slow Assembler
===

A de novo genome assembler written by Peter Anselmo. Written using C, C++, STL and OpenMP.

This assembler using a two-stage greedy assembly algorithm. The first stage trims the reads and assembles into contigs only allowing exact matches. The second stage trims the contigs, and assembles using the Smith Waterman alogrithm. This allows for proper resolution of SNPs and (small) Indels.

This assembler has been developed, tested and optimized for using Illumina HiSeq data and small genomes.

###System Requirements
* A recent version of gcc with support for the OpenMP Library

####Compilation
* Download and extract files, move to directory, and type `make`

###Usage
The program expects exactly one input, a fastq file. a new file "out.fasta" will be generated after the program completes.  Debugging information will be printed to standard out.

```
$ ./ssa sample_reads/phix_100k.fastq > output.txt
```

The parameters for trimming, assembly, and more are located in the `settings.cpp` file. These can be modified to your preference.  Recompilation is required after modification.

###Feedback
Please let me know if you spot anything that can be improved in the code, I am by no means a C/C++ expert.
