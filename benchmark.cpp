#include <iostream>
#include <random>
#include <chrono>
#include <cstring>

#include<fstream>

#include "JPMA_BT.hpp"
#include <time.h>

#define InsertSize 10737418

using namespace std;

void printArguments(){
    cout<<"USAGE: ./benchmark [options]"<<endl;
    cout<<"Options:"<<endl;
    cout<<"    -i [int]     number of key-value pairs to insert"<<endl;
    cout<<"    -d [int]     number of key-value pairs to delete"<<endl;
    cout<<"    -r [int]     length of range for sacnning "<<endl;
    cout<<"    -s [int]     number of key-value pairs to search"<<endl;
    cout<<endl;
}

int main(int argc, char **argv){
    //Redirect cout to file out.txt
    std::ofstream out("out.txt");
    //std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
    
    if (argc == 1) {
        printArguments();
        return 1;
    }

    type_t totalInsert = 0;
    type_t totalDelete = 0;
    type_t rangeLength = 0;
    type_t totalSearch = 0;

    for (type_t i = 1; i<argc; i++) {
        if(strcmp(argv[i], "-i") == 0) {
            totalInsert = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-d") == 0) {
            totalDelete = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-r") == 0) {
            rangeLength = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-s") == 0) {
            totalSearch = atoi(argv[++i]);
        } else {
            printArguments();
            return 1;
        }
    }

    PMA pma;

    if(totalInsert == 0){
        totalInsert = InsertSize;
    }

    if(totalInsert < rangeLength) {
        cout<<"Range length greater than total elements"<<endl;
        exit(1);
    }

    if(totalInsert < totalDelete) {
        cout<<"Number of elements to be deleted greater than total elements"<<endl;
        exit(1);
    }

    int insertCount = 128, inserted = 0;
    int64_t insertDelay = 0;
    chrono::time_point<std::chrono::high_resolution_clock> start, stop;
    while(inserted+insertCount <= totalInsert){
        int64_t *data;
        data =  (int64_t *)malloc(insertCount * sizeof(int64_t));
        srand (time(NULL));
        for(int i = 0; i< insertCount; i++){
            data[i] = i+inserted+1;
            //data[i] = rand() % 1000000;
        }

        for(int i = 0; i<insertCount; i++){
            int64_t source, dest, buffer;
            source = rand() % insertCount;
            dest = rand() % insertCount;
            buffer = data[source];
            data[source] = data[dest];
            data[dest] = buffer;
        }

        start = chrono::high_resolution_clock::now();
        for(int i = 0; i<insertCount; i++){
            pma.insert(data[i], data[i] * 10);
        }
        stop = chrono::high_resolution_clock::now();
        int64_t delay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        if(inserted == 0){
            inserted = insertCount;
        }else{
            inserted += insertCount;
            insertCount *= 2;
        }
        insertDelay += delay;
    }
    cout << "Time taken for insert: " << insertDelay << endl;
    pma.printStat();

    //Searching in the PMA
    type_t records[totalSearch+1];
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> numbers(1,totalSearch);
    for(type_t i=0; i<totalSearch; i++){
        records[i] = numbers(rng);
    }
    start = chrono::high_resolution_clock::now();
    for(type_t i=0; i<totalSearch; i++){
        if(!pma.lookup(records[i])){
            cout<<"Could not get key: "<<records[i]<<endl;
            exit(0);
        }
    }
    stop = chrono::high_resolution_clock::now();
    int64_t searchDelay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    cout<<"Searched "<<totalSearch<<" elements in "<<searchDelay<<" microSeconds."<<endl;
    
    //Range Scan in the PMA
    type_t startRange = rand()%(totalInsert - rangeLength);
    type_t sum_key, sum_value;
    start = chrono::high_resolution_clock::now();
    tie(sum_key, sum_value) = pma.range_sum(startRange, startRange+rangeLength);
    stop = chrono::high_resolution_clock::now();
    if(sum_key*10 != sum_value){
        cout<<"Error in range scan!"<<endl;
        exit(0);
    }
    int64_t scanDelay = chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    cout<<"Scanned elements with range "<<rangeLength<<" in " <<scanDelay<<" microSeconds."<<endl;
    return 0;
}