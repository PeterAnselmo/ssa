#CC = g++ -g -O0 -Wall #debugging
CC = g++ -O3 -Wall
#NVCC = nvcc -O3
NVCC = nvcc -O0 --compiler-options="-g -Wall"
all:ssa
#tests:tests.cpp read.cu sw_matrix.cpp contig.cu fastqfile.cu fasta.cu assembly.cu samfile.cpp
#	$(CC) -o tests tests.cpp
#ssa: ssa.cu read.cu fastqfile.cu fasta.cu contig.cu assembly.cu sw_matrix.cpp
ssa: ssa.cpp settings.cpp read.cpp fastqfile.cpp fasta.cpp contig.cpp assembly.cpp sw_matrix.cpp
	$(CC) -o ssa ssa.cpp
clean:
	rm -f *.o
