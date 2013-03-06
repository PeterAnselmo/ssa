CC = g++ --std=c++0x -O3 -Wall
all:ssa tests
tests:tests.cpp fastqfile.cpp assembly.cpp samfile.cpp
	$(CC) -o tests tests.cpp
ssa:ssa.cpp fastqfile.cpp assembly.cpp samfile.cpp
	$(CC) -o ssa ssa.cpp
