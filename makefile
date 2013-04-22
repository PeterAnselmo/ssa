CC = g++ --std=c++0x -g -O0 -Wall #debugging
#CC = nvcc --std=c++0x -O3 -Wall
NVCC = nvcc -O3
all:ssa
#tests:tests.cpp read.cu sw_matrix.cpp contig.cu fastqfile.cu fasta.cu assembly.cu samfile.cpp
#	$(CC) -o tests tests.cpp
ssa: ssa.cu read.o fastqfile.o fasta.o contig.o assembly.o sw_matrix.cpp
	$(NVCC) -o ssa ssa.cu
read.o:
	$(NVCC) -c read.cu
fastqfile.o:
	$(NVCC) -c fastqfile.cu
fasta.o:
	$(NVCC) -c fasta.cu
contig.o:
	$(NVCC) -c contig.cu
assembly.o:
	$(NVCC) -c assembly.cu
clean:
	rm -f *.o
