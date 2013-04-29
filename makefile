#CC = g++ -g -O0 -Wall #debugging
CC = g++ -O3 -Wall
#NVCC = nvcc -arch sm_30 -O0 --compiler-options="-g -Wall"
NVCC = nvcc -arch sm_30 -O3 --compiler-options="-Wall"
all:ssa
#tests:tests.cu read.cu sw_matrix.cu contig.cu fastqfile.cu fasta.cu assembly.cu samfile.cu
#	$(CC) -o tests tests.cu
#ssa: ssa.cu read.cu fastqfile.cu fasta.cu contig.cu assembly.cu sw_matrix.cu
ssa: ssa.cu settings.cu cudastring.cu read.cu fastqfile.cu fasta.cu contig.cu assembly.cu sw_matrix.cu
	$(NVCC) -o ssa ssa.cu
clean:
	rm -f *.o
