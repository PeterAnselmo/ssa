#CC = g++ --std=c++0x -g -O0 -Wall #debugging
CC = g++ --std=c++0x -O3 -Wall
all:ssa tests
tests:tests.cpp contig.cpp fastqfile.cpp assembly.cpp samfile.cpp
	$(CC) -o tests tests.cpp
ssa: ssa.cpp contig.cpp fastqfile.cpp assembly.cpp samfile.cpp
	$(CC) -I . -o ssa ssa.cpp samtools-0.1.19/libbam.a samtools-0.1.19/bcftools/libbcf.a -lz -lcurses -lpthread 
