CC=g++
CFLAGS=-Wall -g -o3 -std=c++17
INCLUDES=-I ./include/
ALLOC_DEP=./lib/libjemalloc.a
ALLOC_LINK=$(ALLOC_DEP) -lpthread -ldl

PROGRAMS = benchmark

all: $(PROGRAMS)

#btree: 
#	$(CC) $(INCLUDES) $(CFLAGS) -c BPlusTree.cpp -o bptree.o 

jpma:
	$(CC) $(INCLUDES) $(CFLAGS) -c JPMA_BT.cpp -o jpma.o 

benchmark: jpma
	$(CC) $(INCLUDES) $(CFLAGS) jpma.o benchmark.cpp -o benchmark $(ALLOC_LINK)

clean:
	rm -f benchmark jpma.o out.txt