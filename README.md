SSA - Super Slow Assembler
===

A De Novo genome assembler written by Peter Anselmo.  Written using C, C++, STL and OpenMP.  This assembler has been developed, tested and optimized for using Illumina HiSeq data and small genomes.

###System Requirements
* A recent version of gcc with support for the OpenMP Library

####Compilation
* Download and extract files, move to directory, and type `make`

###Usage
The program expects exactly one input, a fastq file. a new file "out.fasta" will be generated after the program completes.  Debugging information will be printed to standard out.

```
$ ./ssa sample_reads/phix_100k.fastq > output.txt
```

###Feedback
Please let me know if you spot anything that can be improved in the code, I am by no means a C/C++ expert.
