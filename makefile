CC = g++ --std=c++0x -O3 -Wall
ssa:ssa.cpp fastafile.cpp assembly.cpp samfile.cpp
	$(CC) -o ssa ssa.cpp
