#include <iostream>
#include <vector>
#include <random>
#include <string.h>
#include <sys/mman.h>
#include <tuple>

#include "defines.hpp"
#include "JPMA_BT.hpp"
//#include "BPlusTree.hpp"
#include "diymalloc.h"

using namespace std;

int treeLevel = 0, leafCount = 0;

PMA::PMA(){
    smallest.push_back(1);                                   //First segment has smallest element 1
    cardinality.push_back(0);                                //First segment contains 0 elements
    totalSegments = 1;                                       //One segment deployed at the start
    lastElementPos.push_back(0);                             //Position of last element in the segment
    elementsInSegment = SEGMENT_SIZE/sizeof(type_t);
    lastValidPos = elementsInSegment - 1;
    MaxThreshold = (elementsInSegment*90)/100;
    MinThreshold = (elementsInSegment*10)/100;
    blocksInSegment = elementsInSegment / JacobsonIndexSize;
    freeSegmentCount = 0;
    type_t *starting_key_chunk, *starting_value_chunk;
    tie(starting_key_chunk, starting_value_chunk) = getSegment();
    
    key_chunks.push_back(starting_key_chunk);
    value_chunks.push_back(starting_value_chunk);
    tree = new BPlusTree(this);
    tree->insertInTree(0, 0, this); //(segment no, dummy key, current JPMA object)
    vector <u_short> blocks;
    for(int i = 0; i<blocksInSegment; i++){
        blocks.push_back(0);
    }
    bitmap.push_back(blocks);

    //Create jacobson Index
    preCalculateJacobson();
}

PMA::~PMA(){
    for(u_int i = 0; i<cleanSegments.size(); i++){
        delete cleanSegments.back();
        cleanSegments.pop_back();
    }
}

void PMA::preCalculateJacobson(){
    for(int i = 0; i<JacobsonIndexCount; i++){
        int j, k = 1, idx = 0;
        u_char count = 0;
        for(j = 0; j<JacobsonIndexSize; j++){
            if(i & k) count++;
            k = k<<1;
        }
        NonZeroEntries[i][idx++] = count;
        k = 1;
        for(j = 0; j<JacobsonIndexSize; j++){
            if(i & k) NonZeroEntries[i][idx++] = j;
            k = k<<1;
        }
    }
}

int PMA::searchSegment(type_t key){
    return tree->searchSegment(key);
}

tuple<type_t *, type_t *> PMA::getSegment(){
    type_t *new_key_chunk, *new_value_chunk;
    
    if(UNLIKELY(freeSegmentCount < 1)){
        if(Allocation_type == 1){
            new_key_chunk = (type_t *) mmap(ADDR, CHUNK_SIZE, PROTECTION, FLAGS, -1, 0);
            if(new_key_chunk == MAP_FAILED){ 
                cout<<"Cannot allocate the virtual memory: " << CHUNK_SIZE << " bytes. mmap error: " << strerror(errno) << "(" << errno << ")"; 
                exit(0);
            }
            new_value_chunk = (type_t *) mmap(ADDR, CHUNK_SIZE, PROTECTION, FLAGS, -1, 0);    
            if(new_value_chunk == MAP_FAILED){ 
                cout<<"Cannot allocate the virtual memory: " << CHUNK_SIZE*16 << " bytes. mmap error: " << strerror(errno) << "(" << errno << ")"; 
                exit(0);
            }
        }
        else{
            new_key_chunk = (type_t *) malloc (CHUNK_SIZE);
            new_value_chunk = (type_t *) malloc (CHUNK_SIZE);    
        }
        cleanSegments.push_back(new_key_chunk);
        cleanSegments.push_back(new_value_chunk);

        freeSegmentCount = CHUNK_SIZE / SEGMENT_SIZE;
        for(int i = 1; i < freeSegmentCount; i++){
            freeKeySegmentBuffer.push_back(new_key_chunk + i * elementsInSegment);
            freeValueSegmentBuffer.push_back(new_value_chunk + i * elementsInSegment);
        }
        freeSegmentCount--;
    }else{
        new_key_chunk = freeKeySegmentBuffer.back();
        new_value_chunk = freeValueSegmentBuffer.back();
        freeKeySegmentBuffer.pop_back();
        freeValueSegmentBuffer.pop_back();
        freeSegmentCount--;
    }
    return {new_key_chunk, new_value_chunk};
}

bool PMA::insert(type_t key, type_t value, int count){
    //Find the location using Binary Search.
    int targetSegment = tree->searchSegment(key);
    
    type_t position = findLocation(key, targetSegment);
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey == key) return false;
    //cout<<"Got location: "<<position<<" Segment: "<<targetSegment<<" cardinality: "<<cardinality[targetSegment]<<" for Key: "<<key<<endl;

    //Check if the current location is empty
    int blockNo = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    if((bitmap[targetSegment][blockNo] & mask) == 0){
        insertInPosition(position, targetSegment, key, value);
        if(cardinality[targetSegment] > (tree->level[0]*SEGMENT_SIZE/8)) tree->redistributeInsert(targetSegment, smallest[targetSegment], this);
        return true;
    }

    //check if need traversing from backside
    if(position >= lastElementPos[targetSegment]){
        if(!insertAfterLast(position, key, value, targetSegment, foundKey, count)) return false;
        if(cardinality[targetSegment] > (tree->level[0]*SEGMENT_SIZE/8)) tree->redistributeInsert(targetSegment, smallest[targetSegment], this);
        return true;
    }

    //Insert among other inserted elements 
    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
    //type_t * segmentKeyOffset = key_chunks[targetSegment];
    int pBase = blockNo * JacobsonIndexSize;
    int insertPos = pBase;
    while(true){
        if(UNLIKELY(ar[0] == 0)){//Got a complete emply block. select the first cell. This can not be the starting block
            insertPos = pBase;
            break;
        }
        else if(ar[0] == JacobsonIndexSize){
            blockNo++;
            if(blockNo == blocksInSegment){
                if(!backSearchInsert(position, key, value, targetSegment, lastValidPos+position)) return false;
                if(cardinality[targetSegment] > (tree->level[0]*SEGMENT_SIZE/8)) tree->redistributeInsert(targetSegment, smallest[targetSegment], this); 
                return true;
            }
            ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
            pBase += JacobsonIndexSize;
            insertPos = pBase;
        }else{
            //Before first element of the block
            if(ar[1] != 0 && position <= pBase){
                insertPos = pBase;
                break;
            } 
            //In between other inserted elements in the block
            bool got = false;
            for(int j = 2; j <= ar[0]; j++){
                if(ar[j]-ar[j-1] > 1 && pBase+ar[j-1] >= position) {
                    insertPos += ar[j-1] + 1;
                    got = true;
                    break;
                }
            }
            //after last element
            if(UNLIKELY(!got)){
                if(ar[ar[0]] == JacobsonIndexSize-1){
                    blockNo++;
                    if(blockNo == blocksInSegment){
                        if(!backSearchInsert(position, key, value, targetSegment, lastValidPos+position)) return false;
                        if(cardinality[targetSegment] > (tree->level[0]*SEGMENT_SIZE/8)) tree->redistributeInsert(targetSegment, smallest[targetSegment], this);
                        return true;
                    }
                    ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
                    pBase += JacobsonIndexSize;
                    insertPos = pBase;
                }else{
                    insertPos += ar[ar[0]] + 1;
                    break;
                }
            }else break;
        }
    }

    if(LIKELY(!backSearchInsert(position, key, value, targetSegment, insertPos))){
        cout<<"error in inserting"<<endl;
        exit(0);
    }
    if(cardinality[targetSegment] > (tree->level[0]*SEGMENT_SIZE/8)) tree->redistributeInsert(targetSegment, smallest[targetSegment], this);
    return true;
}

bool PMA::insertForward(type_t position, type_t key, type_t value, int targetSegment, int insertPos){
    /*
    if(UNLIKELY(insertPos > lastValidPos)) {
        cout<<" No place found for inserting"<<endl;
        exit(0);
    }
    */
    //cout<<"inserting forward. position: "<<position<<" final pos: "<<insertPos<<endl;
    insertInPosition(insertPos, targetSegment, key, value);

    type_t * movePos = key_chunks[targetSegment] + insertPos;
    while(*movePos < *(movePos-1)){
        swapElements(targetSegment, insertPos-1, 1);
        insertPos--;
        if(insertPos == position) break;
        movePos--;
    }
    return true;
}

bool PMA::insertBackward(type_t position, type_t key, type_t value, int targetSegment, int insertPos){
    //cout<<"inserting backward. position: "<<position<<" final pos: "<<insertPos<<endl;
    type_t * segmentKeyOffset = key_chunks[targetSegment];
    insertInPosition(insertPos, targetSegment, key, value);   
    type_t * moveback = segmentKeyOffset + insertPos;
    while(*moveback > *(moveback+1)){
        swapElements(targetSegment, insertPos, 1);
        insertPos++;
        if(insertPos == position) break;
        moveback++;
    }
    return true;
}

bool PMA::insertAfterLast(type_t position, type_t key, type_t value, int targetSegment, type_t foundKey, int count) {
    //type_t * segmentKeyOffset = key_chunks[targetSegment];
    if(UNLIKELY(position == lastValidPos)){ //Got out of the current segment. Traverse backward for vacant space
        //return backSearchInsert(position, key, value, targetSegment, lastValidPos+position);
        for(type_t i = lastValidPos-1; i>=0; i--){
            int blockNo = i/JacobsonIndexSize;
            int bitPosition = i % JacobsonIndexSize;
            u_short mask =  1 << bitPosition;

            if((bitmap[targetSegment][blockNo] & mask) == 0) {
                return insertBackward(position, key, value, targetSegment, i);
            }
        }
        cout<<"Program should never reach hear at InsertAfterLast"<<endl;
        exit(0);
    }
    else{ //Have some space left in the segment. Go forward in the space max 3 slots
        type_t adjust = min(lastValidPos - lastElementPos[targetSegment], (type_t) MaxGap);
        adjust = min(adjust, abs(*(key_chunks[targetSegment]+lastElementPos[targetSegment])-key));
        if(key > foundKey){
            //cout<<"inserting forward. position: "<<position<<" final pos: "<<position+adjust<<endl;
            insertInPosition(position+adjust, targetSegment, key, value);
            return true;
        }
        else{
            //cout<<"inserting forward. position: "<<position<<" final pos: "<<position+adjust<<endl;
            insertInPosition(position+adjust, targetSegment, key, value);
            swapElements(targetSegment, position, adjust);
            return true;
        }
    }
}

bool PMA::backSearchInsert(type_t position, type_t key, type_t value, int targetSegment, int forwardInsertPos) {
    type_t blockNo = position/JacobsonIndexSize; 
    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
    type_t insertPos = blockNo * JacobsonIndexSize;
    
    while(true){
        if(ar[0] == 0){//Got a complete empty block. select the last cell
            insertPos += JacobsonIndexSize - 1;
            break;
        }
        else if(ar[0] == JacobsonIndexSize){ //Fully loaded block
            blockNo--;
            if(blockNo < 0){
                return insertForward(position, key, value, targetSegment, forwardInsertPos);
            }
            ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
            insertPos -= JacobsonIndexSize;
        }else{
            int elements = ar[0];
            bool got = false;
            if(ar[elements] < (JacobsonIndexSize - 1) && insertPos + ar[elements] < position){
                insertPos += JacobsonIndexSize - 1;
                got = true;
            }else{
                for(int j = elements; j > 1; j--){
                    if(ar[j]-ar[j-1] > 1 && insertPos + ar[j] <= position) {
                        insertPos += ar[j] - 1;
                        got = true;
                        break;
                    }
                }
            }
            if(!got){
                if(ar[1] == 0){
                    blockNo--;
                    if(blockNo < 0){
                        return insertForward(position, key, value, targetSegment, forwardInsertPos);
                    }
                    
                    ar = NonZeroEntries[bitmap[targetSegment][blockNo]];
                    insertPos -= JacobsonIndexSize;
                }else{
                    insertPos += ar[1] - 1;
                    break;
                }
            }else break;
        }

    }
    //if(insertPos < 0) insertPos += lastValidPos;
    if(abs(insertPos-position)>abs(forwardInsertPos-position)){
        return insertForward(position, key, value, targetSegment, forwardInsertPos);
    }else return insertBackward(position, key, value, targetSegment, insertPos);
}

void PMA::swapElements(type_t targetSegment, type_t position, type_t adjust){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t holdKey = *(segmentOffset + position);
    *(segmentOffset + position) = *(segmentOffset + position + adjust);
    *(segmentOffset + position + adjust) = holdKey;

    segmentOffset = value_chunks[targetSegment];
    type_t holdValue = *(segmentOffset + position);
    *(segmentOffset + position) = *(segmentOffset + position + adjust);
    *(segmentOffset + position + adjust) = holdValue;
}

void PMA::insertInPosition(type_t position, int targetSegment, type_t key, type_t value){
    //Store key, value and update bitmap, cardinality and last index
    type_t * segmentOffset = key_chunks[targetSegment];
    *(segmentOffset + position) = key;
    segmentOffset = value_chunks[targetSegment];
    *(segmentOffset + position) = value;
    int blockPosition = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    bitmap[targetSegment][blockPosition] |= mask;
    cardinality[targetSegment]++;
    if(lastElementPos[targetSegment] < position) lastElementPos[targetSegment] = position;
    //if(cardinality[targetSegment] > MaxThreshold) redistribute(targetSegment, INSERT);
}

bool PMA::remove(type_t key){
    int targetSegment = searchSegment(key);

    type_t position = findLocation(key, targetSegment);
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t foundKey = *(segmentOffset + position);
    if(foundKey != key) return false;
    deleteInPosition(position, targetSegment, key);
    return true;
}

void PMA::deleteInPosition(type_t position, int targetSegment, type_t key){
    int blockPosition = position/JacobsonIndexSize;
    int bitPosition = position % JacobsonIndexSize;
    u_short mask =  1 << bitPosition;
    bitmap[targetSegment][blockPosition] &= (~mask);
    cardinality[targetSegment]--;
    if(lastElementPos[targetSegment] == position){
        u_char * ar = NonZeroEntries[bitmap[targetSegment][blockPosition]];
        lastElementPos[targetSegment] = blockPosition * JacobsonIndexSize + ar[ar[0]];
    }

    //Will be handled later
    //if(cardinality[targetSegment] < MinThreshold) redistribute(targetSegment, DELETE);
}

void PMA::deleteSegment(int targetSegment){
    freeSegmentCount++;
    freeKeySegmentBuffer.push_back(key_chunks[targetSegment]);
    freeValueSegmentBuffer.push_back(value_chunks[targetSegment]);

    key_chunks.erase(key_chunks.begin() + targetSegment);
    value_chunks.erase(value_chunks.begin() + targetSegment);
    smallest.erase(smallest.begin() + targetSegment);
    lastElementPos.erase(lastElementPos.begin() + targetSegment);
    cardinality.erase(cardinality.begin() + targetSegment);
    bitmap.erase(bitmap.begin() + targetSegment);
    totalSegments--;
}

bool PMA::lookup(type_t key){
    int targetSegment = tree->searchSegment(key);

    type_t position = findLocation(key, targetSegment);
    type_t * segmentOffsetKey = key_chunks[targetSegment];
    type_t * segmentOffsetVal = value_chunks[targetSegment];
    type_t foundKey = *(segmentOffsetKey + position);
    type_t foundVal = *(segmentOffsetVal + position);
    if(foundKey*10 != foundVal) {cout<<"error in the tree while searching"<<endl; exit(0);}
    if(foundKey == key) return true;  
    return false;
}

type_t PMA::findLocation(type_t key, int targetSegment){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t start = 0, end = lastElementPos[targetSegment];
    int blockPosition, bitPosition, mask;
    type_t data, mid = 0;
    while(start <= end){
        mid = (start + end) / 2;
        blockPosition = mid / JacobsonIndexSize;
        bitPosition = mid % JacobsonIndexSize;
        mask = 1 << bitPosition;
        if((bitmap[targetSegment][blockPosition] & mask) == 0){
            int64_t changedMid = mid, offset = -1;
            while(changedMid >= start){
                changedMid += offset;
                blockPosition = changedMid / JacobsonIndexSize;
                bitPosition = changedMid % JacobsonIndexSize;
                mask = 1 << bitPosition;
                if((bitmap[targetSegment][blockPosition] & mask) !=0) break;
            }
            if(changedMid < start){
                changedMid = mid;
                offset = 1;
                while(changedMid <= end){
                    changedMid += offset;
                    blockPosition = changedMid / JacobsonIndexSize;
                    bitPosition = changedMid % JacobsonIndexSize;
                    mask = 1 << bitPosition;
                    if((bitmap[targetSegment][blockPosition] & mask) !=0 ) break;
                }
                if(changedMid > end){
                    return mid;
                }else mid = changedMid;
            }else mid = changedMid;
        }
        data = *(segmentOffset + mid);
        if(data == key) return mid;
        else if(data < key) start = mid + 1;
        else end = mid - 1;
    }
    return mid;
}

type_t PMA::findLocation2(type_t key, int targetSegment){
    type_t * segmentOffset = key_chunks[targetSegment];
    type_t start = 0;
    type_t end = lastElementPos[targetSegment];
    type_t data, mid = 0;
    int blockPosition, bitPosition, mask;
    while(start<=end){
        mid = (start+end)/2;
        blockPosition = mid / JacobsonIndexSize;
        bitPosition = mid % JacobsonIndexSize;
        mask = 1 << bitPosition;
        if((bitmap[targetSegment][blockPosition] & mask) == 0){
            type_t base = blockPosition * JacobsonIndexSize;
            bool over = false, found = false;
            while(true){
                u_char * ar = NonZeroEntries[bitmap[targetSegment][blockPosition]];
                if(UNLIKELY(ar[0] == 0)){
                    blockPosition++;
                    base += JacobsonIndexSize;
                    if(blockPosition >= blocksInSegment || base > end) break;
                    continue;
                }
                for(int i = 1; i<=ar[0]; i++){
                    if(base + ar[i] > mid) {
                        if(base + ar[i] <= end) {mid = base + ar[i]; found = true;}
                        else over = true;
                        break;
                    }
                }
                if(found || over) break;
                blockPosition++;
                base += JacobsonIndexSize;
            }
            if(!found){
                blockPosition = mid / JacobsonIndexSize;
                type_t base = blockPosition * JacobsonIndexSize;
                while(true){
                    u_char * ar = NonZeroEntries[bitmap[targetSegment][blockPosition]];
                    if(UNLIKELY(ar[0] == 0)){
                        blockPosition--;
                        base -= JacobsonIndexSize;
                        if(blockPosition < 0 || base < start) return mid;
                        continue;
                    }
                    for(int i = ar[0]; i>0; i--){
                        if(base + ar[i] < mid){
                            if(base + ar[i] < start) return mid;
                            mid = base + ar[i];
                            found = true;
                            break;
                        }
                    }
                    if(found) break;
                    blockPosition--;
                    base -= JacobsonIndexSize;
                }
            }
        }
        data = *(segmentOffset + mid);
        if(data == key) return mid;
        else if(data < key) start = mid + 1;
        else end = mid - 1;
    }
    return mid;
}

void PMA::printAllElements(){
    tree->printAllElements(this);
}

void PMA::printSegElements(int targetSegment){
    type_t * key = key_chunks[targetSegment];
    type_t pBase = 0;
    for(type_t block = 0; block<blocksInSegment; block++){
        u_short bitpos = 1;
        cout <<" Bitmap: "<<bitmap[targetSegment][block]<<" ";
        for(type_t j = 0; j<JacobsonIndexSize; j++){
            if(bitmap[targetSegment][block] & bitpos)
                cout << *(key+pBase+j) << " ";
            else cout <<"0 ";
            bitpos = bitpos << 1;
        }
        pBase += JacobsonIndexSize;
    }
    cout<<"last offset: "<<lastElementPos[targetSegment] <<" Cardinality: "<<cardinality[targetSegment]<<" Total Segment: "<<totalSegments<< endl;
}

void PMA::printStat(){
    type_t totalElements = 0;
    for(type_t i = 0; i<totalSegments; i++){
        totalElements += cardinality[i];
    }
    cout<<"Total elements: "<<totalElements<<endl;
    cout<<"Total Segment: "<<totalSegments<<", Free Segments: "<<freeSegmentCount<<", Elements in a Segment: "<<elementsInSegment<<endl;
    cout<<"Redistribute with insert: "<<redisInsCount<<", Redistribute with update: "<<redisUpCount<<endl;
    /*
    vector<BPlusTree::node *> temp;
    temp.push_back(tree->root);
    treeLevel = leafCount = 0;
    tree->printTree(temp, 1);
    */
}


BPlusTree::BPlusTree(PMA *obj){
    root = NULL;
    calculateThreshold();
}

void BPlusTree::insertInTree(int chunkNo, type_t search_key, PMA *obj){
    if(UNLIKELY(root == NULL)){
        root = new Node();
        leaf *leafNode = new leaf();
        root->child_ptr[0] = (node *)leafNode;
        root->key[0] = INT64_MAX;
        root->ptrCount = 1;
        root->nodeLeaf = true;

        leafNode->segNo[0] = chunkNo;
        leafNode->key[0] = INT64_MAX;
        leafNode->childCount = 1;
        return;
    }

    leaf *leaf = findLeaf(search_key);
    if(leaf->childCount < Leaf_Degree){ //insert in leaf
        if(leaf->childCount == 1){
            int segNo = leaf->segNo[0];
            if(obj->smallest[segNo] > search_key){
                leaf->segNo[1] = leaf->segNo[0];
                leaf->segNo[0] = chunkNo;
                leaf->key[0] = obj->smallest[segNo];
            }else{
                leaf->segNo[1] = chunkNo;
                leaf->key[0] = search_key;
            }
            leaf->childCount = 2;
            return;
        }
        int position;
        for(position = leaf->childCount-1; position>0; position--){
            if(leaf->key[position-1] > search_key){
                leaf->key[position] = leaf->key[position-1];
                leaf->segNo[position+1] = leaf->segNo[position]; 
            }
            else{
                int segNo = leaf->segNo[position];
                if(obj->smallest[segNo] > search_key){
                    leaf->segNo[position+1] = leaf->segNo[position];
                    leaf->segNo[position] = chunkNo;
                    leaf->key[position] = obj->smallest[segNo];
                }else{
                    leaf->segNo[position+1] = chunkNo;
                    leaf->key[position] = search_key;
                }
                leaf->childCount++;
                /*
                if(leaf->childCount > Leaf_Degree){
                    cout<<"Got more child than limit: 1"<<endl;
                    cout<<"Leaf childcount: "<<leaf->childCount<<endl;
                    exit(0);
                }
                */
                return;
            }
        }
        int segNo = leaf->segNo[0];
        if(obj->smallest[segNo] > search_key){
            leaf->segNo[1] = leaf->segNo[0];
            leaf->segNo[0] = chunkNo;
            leaf->key[0] = obj->smallest[segNo];
        }else{
            leaf->segNo[1] = chunkNo;
            leaf->key[0] = search_key;
        }
        leaf->childCount++;
        /*
        if(leaf->childCount > Leaf_Degree){
            cout<<"Got more child than limit: 2"<<endl;
            cout<<"Leaf childcount: "<<leaf->childCount<<endl;
            exit(0);
        }
        */
        return;
    }
    //if(leaf->childCount < Leaf_Degree){cout<<"Progrom should not have reached here: InsertInTree"<<endl; exit(0);}
    //Leaf is not empty. Divide.
    //Copy the key-value pairs
    type_t key_store[Leaf_Degree+1];
    int segNo_store[Leaf_Degree+1];
    int position;
    bool done = false;
    for(position = leaf->childCount-1; position>0; position--){
        if(leaf->key[position-1] > search_key){
            key_store[position] = leaf->key[position-1];
            segNo_store[position + 1] = leaf->segNo[position];
        }else{
            int segNo = leaf->segNo[position];
            if(obj->smallest[segNo] > search_key){
                segNo_store[position+1] = leaf->segNo[position];
                segNo_store[position] = chunkNo;
                key_store[position] = obj->smallest[segNo];
            }else{
                segNo_store[position+1] = chunkNo;
                segNo_store[position] = leaf->segNo[position];
                key_store[position] = search_key;
            }
            done = true;
            position--;
            break;
        }
    }
    if(done){
        for( ; position>=0; position--){
            segNo_store[position] = leaf->segNo[position];
            key_store[position] = leaf->key[position];
        }
    }else{
        int segNo = leaf->segNo[0];
        if(obj->smallest[segNo] > search_key){
            segNo_store[position+1] = leaf->segNo[position];
            segNo_store[position] = chunkNo;
            key_store[position] = obj->smallest[segNo];
        }else{
            segNo_store[position+1] = chunkNo;
            segNo_store[position] = leaf->segNo[position];
            key_store[position] = search_key;
        }
    }

    //Copying done. Create two nodes
    BPlusTree::leaf *l2 = new BPlusTree::leaf();
    for(int halfLeaf = 0; halfLeaf <= Leaf_Degree/2; halfLeaf++){
        leaf->segNo[halfLeaf] = segNo_store[halfLeaf];
        leaf->key[halfLeaf] = key_store[halfLeaf];
    }
    for(int nextHalf = Leaf_Degree/2 + 1, index = 0; nextHalf <= Leaf_Degree; nextHalf++, index++){
        l2->segNo[index] = segNo_store[nextHalf];
        l2->key[index] = key_store[nextHalf];
    }
    leaf->childCount = Leaf_Degree/2 + 1;
    l2->childCount = Leaf_Degree + 1 - leaf->childCount;
    /*
    if(leaf->childCount > Leaf_Degree){
        cout<<"Got more child than limit: 3"<<endl;
        cout<<"Leaf childcount: "<<leaf->childCount<<endl;
        exit(0);
    }
    if(l2->childCount > Leaf_Degree){
        cout<<"Got more child than limit: 4"<<endl;
        cout<<"Leaf childcount: "<<l2->childCount<<endl;
        exit(0);
    }
    */
    l2->nextLeaf = leaf->nextLeaf;
    leaf->nextLeaf = l2;
    type_t key_parent = obj->smallest[leaf->segNo[0]];
    insert_in_parent(leaf,key_store[Leaf_Degree/2],l2, key_parent);
}

void BPlusTree::insert_in_parent(void *left, type_t search_key, void *right, type_t key_parent){
    if(left == root || right == root){
        node *N = new node();
        N->child_ptr[0] = (node *)left;
        N->child_ptr[1] = (node *)right;
        N->key[0] = search_key; 
        N->ptrCount = 2;
        root = N;
        return;
    }
    node *N = findParent(left, key_parent);
    if(N->ptrCount < Tree_Degree){
        for(int position = N->ptrCount-1; position >= 0; position--){
            if(N->child_ptr[position] == left){
                N->child_ptr[position+1] = (node *)right;
                N->key[position] = search_key;
                N->ptrCount++;
                return;
            }else{
                N->child_ptr[position+1] = N->child_ptr[position];
                N->key[position] = N->key[position-1];
            }
        }
        cout<<"Program should never reach here. Insert in Parent node of B+ Tree"<<endl;
        exit(0);
    }
    /*
    if(N->ptrCount>Tree_Degree){
        cout<<"Non Leaf node has more children than tree degree"<<endl;
        vector<node *> temp;
        temp.push_back(root);
        treeLevel = leafCount = 0;
        printTree(temp, 1);
        exit(0);
    }
    */
    type_t key_store[Tree_Degree+1];
    Node *ptr_store[Tree_Degree+1];
    int position;
    bool done = false;
    for(position = N->ptrCount-1; position >= 0; position--){
        if(N->child_ptr[position] == left){
            ptr_store[position+1] = (node *) right;
            ptr_store[position] = N->child_ptr[position];
            key_store[position] = search_key;
            done = true;
        }else if(done){
            ptr_store[position] = N->child_ptr[position];
            key_store[position] = N->key[position];
        }else{
            ptr_store[position+1] = N->child_ptr[position];
            key_store[position] = N->key[position-1];
        }
    }

    //Spilit records into two nodes
    node *N2 = new node();
    for(int half = 0; half <= Tree_Degree/2; half++){
        N->child_ptr[half] = ptr_store[half];
        N->key[half] = key_store[half];
    }
    for(int nextHalf = Tree_Degree/2 + 1, index = 0; nextHalf <= Tree_Degree; nextHalf++, index++){
        N2->child_ptr[index] = ptr_store[nextHalf];
        N2->key[index] = key_store[nextHalf];
    }
    N->ptrCount = Tree_Degree/2 + 1;
    N2->ptrCount = Tree_Degree + 1 - N->ptrCount;
    N2->nodeLeaf = N->nodeLeaf;
    insert_in_parent(N, key_store[Tree_Degree/2], N2, key_parent);
}

BPlusTree::node* BPlusTree::findParent(void *n, type_t search_key){
    if(n == root) return NULL;
    node *parent = root;
    while(true){
        bool found = false;
        for(int position = parent->ptrCount-1; position>0; position--){
            if(parent->key[position-1] <= search_key){
                if(n == parent->child_ptr[position]) return parent;
                else{
                    parent = parent->child_ptr[position];
                    found = true; 
                    break;
                }
            }
        }
        if(!found){
            if(n == parent->child_ptr[0]) return parent;
            parent = parent->child_ptr[0];
        }
    }
}

BPlusTree::leaf* BPlusTree::findLeaf(type_t search_key){
    node *temp = root;
    while(!temp->nodeLeaf){
        int smallest;
        for(smallest = temp->ptrCount - 1; smallest > 0; smallest--){
            if(temp->key[smallest-1] <= search_key) break;
        }
        temp = temp->child_ptr[smallest];
    }
    for(int smallest = temp->ptrCount - 1; smallest > 0; smallest--){
        if(temp->key[smallest-1] <= search_key) return (leaf *)temp->child_ptr[smallest];
    }
    return (leaf *)temp->child_ptr[0];
}

int BPlusTree::searchSegment(type_t search_key){
    leaf *leaf = findLeaf(search_key);
    if(leaf->childCount == 1) return leaf->segNo[0];
    for(int smallest = leaf->childCount - 1; smallest > 0; smallest--){
        if(leaf->key[smallest-1] <= search_key) return leaf->segNo[smallest];
    }
    return leaf->segNo[0];
}

void BPlusTree::calculateThreshold(){
    level[0] = 0.95;
    level[MaxLevel] = 0.50;
    for(int i = 1; i<MaxLevel; i++){
        level[i] = level[0] - ((level[0]-level[MaxLevel])/(MaxLevel-i));
    }
}

void BPlusTree::redistributeInsert(int segment, type_t SKey, PMA *obj){
    obj->redisInsCount++;
    leaf *par = findLeaf(SKey);
    if(findCardinality(par, obj) >= (level[1]*Leaf_Degree*SEGMENT_SIZE/8)){
        cout<<"Got upper level of tree"<<endl;
        node *parent = findParent(par, SKey);
        type_t nodeCard = findCardinality(parent, obj);
        int tree_Degree_Count = Leaf_Degree * Tree_Degree;
        if(nodeCard >= level[2]*tree_Degree_Count*SEGMENT_SIZE/8){
            int cLevel = 2;
            while(nodeCard >= level[cLevel]*tree_Degree_Count*SEGMENT_SIZE/8){
                cLevel++;
                tree_Degree_Count *= Tree_Degree;
                parent = findParent(parent, SKey);
                nodeCard = findCardinality(parent, obj);
                if (parent == root) break;
            }
            vector<int> segments;
            listSegments(segments, parent);
            //Logically delete the branch
            leaf *left = leftmostLeaf(parent);
            left = findLeaf(obj->smallest[left->segNo[0]]);
            leaf *right = rightmostLeaf(parent);
            left->nextLeaf = right->nextLeaf;
            //Physically delete
            deleteNode(parent);
    
            reinsertInTree(segments, nodeCard, obj);
        }else{ //Only redistribute the segments under the leaf node
            nodeCard = findCardinality(par, obj);
            vector<int> segments;
            for(int i=0; i < par->childCount; i++){
                segments.push_back(par->segNo[i]);
            }
            
            if(obj->smallest[par->segNo[0]] != 0){
                deleteLeaf(par, obj->smallest[par->segNo[0]]);
                //delete par;
            }
            reinsertInTree(segments, nodeCard, obj);
        }
    }else{
        //Divide in 2 segments
        int segNo = obj->redistributeWithDividing(segment);
        insertInTree(segNo, obj->smallest[segNo], obj);
    }
}

void BPlusTree::reinsertInTree(vector<int> &segments, type_t cardi, PMA *obj){
    vector<type_t *> usedKeySegments;
    vector<type_t *> usedValSegment;
    vector<type_t> smallest;
    vector<type_t> lastElementPos;
    vector<int> cardinality;
    vector<vector<u_short>> bitmap;

    type_t insertPos = 0;
    type_t prevKey;
    type_t lowestElement = 0;
    int elementCount = 0;
    type_t *keyStore, *valueStore;
    bool isEmptySegment = true;
    type_t maxGap = (type_t) MaxGap;
    vector <u_short> blocks;

    tie(keyStore, valueStore) = obj->getSegment();
    for(int i = 0; i<obj->blocksInSegment; i++){
        blocks.push_back(0);
    }

    cout<<"Inside redistribution. Segments: "<<segments.size()<<" Cardinality: "<<cardi<<" Segments are: ";
    for(u_int i = 0; i < segments.size(); i++){
        cout<<segments[i]<<" ";
    }cout<< " Segment Limit: "<<level[0]*SEGMENT_SIZE/128<<endl;

    //Copy the content to new set of segments
    for(u_int i = 0; i<segments.size(); i++){
        int segNo = segments[i];
        type_t *sourceKey = obj->key_chunks[segNo];
        type_t *sourceVal = obj->value_chunks[segNo];
        type_t position = 0;
        for(int bl = 0; bl <obj->blocksInSegment; bl++ ){
            u_char * ar = obj->NonZeroEntries[obj->bitmap[segNo][bl]];
            for(int j = 1; j <= ar[0]; j++){
                if(isEmptySegment) {
                    *keyStore = prevKey = *(sourceKey + ar[j]);
                    *valueStore = *(sourceVal + ar[j]);
                    lowestElement = prevKey;
                    elementCount++;
                    isEmptySegment = false;
                }else{
                    type_t curKey = *(sourceKey + ar[j]);
                    int gap = min(curKey - prevKey, maxGap);
                    prevKey = curKey;
                    if(gap < 0) {
                        cout<<endl<<"Curlock " <<position + ar[j]<<endl;
                        cout<<"curKey: "<<curKey<<" PrevKey: "<<prevKey<<endl;
                        cout<<"Gap should not be negative"<<endl;
                        printAllElements(obj);
                        exit(0);
                        //*keyStore = *(sourceKey-100000000);
                    }
                    if(insertPos + gap > obj->lastValidPos || elementCount == level[0]*SEGMENT_SIZE/8){
                        //save prev seg
                        usedKeySegments.push_back(keyStore);
                        usedValSegment.push_back(valueStore);
                        smallest.push_back(lowestElement);
                        lastElementPos.push_back(insertPos);
                        cardinality.push_back(elementCount);
                        bitmap.push_back(blocks);

                        //create new seg and reset variables
                        insertPos = 0;
                        lowestElement = 0;
                        elementCount = 0;
                        isEmptySegment = true;
                        tie(keyStore, valueStore) = obj->getSegment();
                        vector <u_short> newBlock;
                        for(int ii = 0; ii<obj->blocksInSegment; ii++){
                            newBlock.push_back(0);
                        }
                        
                        cout<<"Before changing"<<endl;
                        for(int ii = 0; ii<obj->blocksInSegment; ii++){
                            cout<<" "<<blocks[ii];
                        }
                        blocks = newBlock;
                        for(int ii = 0; ii<obj->blocksInSegment; ii++){
                            cout<<" "<<blocks[ii];
                        }
                        j--;
                    }else{
                        insertPos += gap;
                        *(keyStore + insertPos) = curKey;
                        *(valueStore + insertPos) = *(sourceVal + ar[j]);
                        int blockPosition = insertPos / JacobsonIndexSize;
                        int bitPosition = insertPos % JacobsonIndexSize;
                        u_short mask = 1 << bitPosition;
                        blocks[blockPosition] |= mask;
                        elementCount++;
                    }
                }
            }
            sourceKey += JacobsonIndexSize;
            sourceVal += JacobsonIndexSize;
            position += JacobsonIndexSize;
        }
    }

    //place the created segments to the final group of segments
    for(u_int i = 0; i < segments.size(); i++){
        obj->freeKeySegmentBuffer.push_back(obj->key_chunks[segments[i]]);
        obj->freeValueSegmentBuffer.push_back(obj->value_chunks[segments[i]]);
        obj->freeSegmentCount++;
    }

    u_int j = 0;
    int lastSeg = obj->totalSegments;
    for(u_int i=0; i<usedKeySegments.size(); i++){

        if(j >= segments.size()){
            obj->key_chunks.push_back(usedKeySegments[i]);
            obj->value_chunks.push_back(usedValSegment[i]);
            obj->smallest.push_back(smallest[i]);
            obj->lastElementPos.push_back(lastElementPos[i]);
            obj->cardinality.push_back(cardinality[i]);
            obj->bitmap.push_back(bitmap[i]);
            segments.push_back(lastSeg++);
            obj->totalSegments++;
        }else{
            obj->key_chunks[segments[j]] = usedKeySegments[i];
            obj->value_chunks[segments[j]] = usedValSegment[i];
            obj->smallest[j] = smallest[i];
            obj->lastElementPos[j] = lastElementPos[i];
            obj->cardinality[j] = cardinality[i];
            obj->bitmap[j] = bitmap[i];
            j++;
        }
    }

    //Check all the segments.
    for(u_int i = 0; i <obj->key_chunks.size(); i++)
    {
        cout<<"Key Chunk: "<<(ulong)obj->key_chunks[i]<<endl;
        cout<<"Smallest: "<<obj->smallest[i]<<endl;
        cout<<"LastElementPos: "<<obj->lastElementPos[i]<<endl;
        cout<<"Cardinality: "<<obj->cardinality[i]<<endl;
        cout<<"Bitmap: ";
        for(int j=0; j<obj->blocksInSegment; j++){
            cout<<" "<<obj->bitmap[i][j];
        }
    }

    //Last part. Insert the segments to the tree
    for(u_int i=0; i<segments.size(); i++){
        insertInTree(segments[i], obj->smallest[segments[i]], obj);
    }
}

BPlusTree::leaf* BPlusTree::rightmostLeaf(node *parent){
    if(parent->nodeLeaf) return (leaf *)parent->child_ptr[parent->ptrCount-1];
    while(!parent->nodeLeaf){
        parent = parent->child_ptr[root->ptrCount-1];
    }
    return (leaf *)parent->child_ptr[parent->ptrCount-1];
}

BPlusTree::leaf* BPlusTree::leftmostLeaf(node *parent){
    if(parent->nodeLeaf) return (leaf *)root->child_ptr[0];
    while(!parent->nodeLeaf){
        parent = parent->child_ptr[0];
    }
    return (leaf *)parent->child_ptr[0];
}

void BPlusTree::listSegments(vector<int> &segments, node *parent){
    if(parent->nodeLeaf){
        for(int i = 0; i<parent->ptrCount; i++){
            leaf *l = (leaf *)parent->child_ptr[i];
            for(int j = 0; j<l->childCount; j++){
                segments.push_back(l->segNo[j]);
            }
            delete l;
        }
        return;
    }
    for(int i = 0; i<parent->ptrCount; i++){
        listSegments(segments, parent->child_ptr[i]);
        delete parent->child_ptr[i];
    }
}

void BPlusTree::deleteNode(node *parent){
    if(parent->nodeLeaf){
        for(int i=0; i<parent->ptrCount; i++){
            delete (leaf *)parent->child_ptr[i];
        }
        delete parent;
    }else{
        for(int i = 0; i<parent->ptrCount; i++){
            deleteNode(parent->child_ptr[i]);
            delete parent;
        }
    }
}

void BPlusTree::deleteLeaf(leaf *l, type_t SKey){
    leaf *prev = findLeaf(SKey);
    prev->nextLeaf = l->nextLeaf;

    node *par = findParent(l, SKey);
    while(par->ptrCount == 1){
        if(par == root) return;
        node *par2 = findParent(par, SKey);
        par2->child_ptr[0] = par->child_ptr[0];
        delete par;
        par = par2;
    }

}

/*
    Returns new segment nubmer. Unsed in cases only one new segment needs to be created
 */
int PMA:: redistributeWithDividing(int targetSegment){
    type_t halfElement = cardinality[targetSegment]/2;
    type_t *new_key_chunk, *new_value_chunk;
    tie(new_key_chunk, new_value_chunk) = getSegment();

    type_t * moveKeyOffset = key_chunks[targetSegment];
    type_t * moveValOffset = value_chunks[targetSegment];
    type_t * destKeyOffset = new_key_chunk;
    type_t * destValOffset = new_value_chunk;

    type_t copyBlock, i, j, elementCount = 0;

    for(copyBlock=0; copyBlock<blocksInSegment; copyBlock++){
        u_char * ar = NonZeroEntries[bitmap[targetSegment][copyBlock]];
        if((elementCount + ar[0]) >= halfElement){
            halfElement = elementCount + ar[0];
            lastElementPos[targetSegment] = copyBlock * JacobsonIndexSize + ar[ar[0]];
            copyBlock++;
            break;
        }
        elementCount += ar[0];
    }        
    vector<u_short> blocks;
    
   for(int ii = 0; ii < blocksInSegment; ii++){
        blocks.push_back(0);
    }

    //Copy the elements of current block
    elementCount = cardinality[targetSegment]-halfElement;
    u_char * ar = NonZeroEntries[bitmap[targetSegment][copyBlock]];
    while(UNLIKELY(ar[0] == 0)){
        copyBlock++;
        ar = NonZeroEntries[bitmap[targetSegment][copyBlock]];
    }
    type_t * pKeyBase = moveKeyOffset + copyBlock * JacobsonIndexSize;
    type_t * pValBase = moveValOffset + copyBlock * JacobsonIndexSize;
    type_t lastInput = 0, lastInsertkey = 0;

    *destKeyOffset = lastInsertkey = *(pKeyBase + ar[1]);
    *destValOffset = *(pValBase + ar[1]);
    blocks[0] = 1;
    for(i = 2, j = 0; i<=ar[0]; i++){
        type_t current_element = *(pKeyBase + ar[i]);
        int keyGap = current_element - lastInsertkey;
        if(keyGap < MaxGap) j += keyGap;
        else j+= MaxGap;
        *(destKeyOffset + j) = lastInsertkey = current_element;
        *(destValOffset + j) = *(pValBase + ar[i]);
        int blockPosition = j / JacobsonIndexSize;
        int bitPosition = j % JacobsonIndexSize;
        u_short mask = 1 << bitPosition;
        blocks[blockPosition] |= mask;
    }
    lastInput = j;
    elementCount -= ar[0];
    bitmap[targetSegment][copyBlock] = 0;

    type_t lastAccessPos = lastValidPos - 1;
    for(type_t blockno = copyBlock+1; blockno < blocksInSegment; blockno++){
        pKeyBase += JacobsonIndexSize;
        pValBase += JacobsonIndexSize;
        ar = NonZeroEntries[bitmap[targetSegment][blockno]];
        if(UNLIKELY(ar[0] == 0)){bitmap[targetSegment][blockno] = 0; continue;}
        for(i = 1; i<=ar[0]; i++){
            type_t current_element = *(pKeyBase + ar[i]);
            if(elementCount < (lastAccessPos - j)){
                type_t keyGap = current_element - lastInsertkey;
                if(keyGap<MaxGap) j += keyGap;
                else j+= MaxGap;
            }else j++;
            *(destKeyOffset + j) = lastInsertkey = current_element;
            *(destValOffset + j) = *(pValBase + ar[i]);
            int blockPosition = j / JacobsonIndexSize;
            int bitPosition = j % JacobsonIndexSize;
            u_short mask = 1 << bitPosition;
            blocks[blockPosition] |= mask;
            elementCount--;
        }
        lastInput = j;
        bitmap[targetSegment][blockno] = 0;
    }

    //key_chunks.emplace(key_chunks.begin() + targetSegment + 1, new_key_chunk);
    key_chunks.push_back(new_key_chunk);
    //value_chunks.emplace(value_chunks.begin() + targetSegment + 1, new_value_chunk);
    value_chunks.push_back(new_value_chunk);
    //lastElementPos.emplace(lastElementPos.begin() + targetSegment + 1, lastInput);
    lastElementPos.push_back(lastInput);
    //lastElementPos[targetSegment] = lastElementTargetSegment;
    //smallest.emplace(smallest.begin() + targetSegment + 1, *(destKeyOffset));
    smallest.push_back(*(destKeyOffset));
    //cardinality.emplace(cardinality.begin() + targetSegment + 1, cardinality[targetSegment] - halfElement);
    cardinality.push_back(cardinality[targetSegment] - halfElement);
    cardinality[targetSegment] = halfElement;
    //bitmap.emplace(bitmap.begin() + targetSegment + 1, blocks);
    bitmap.push_back(blocks);
    totalSegments++;
    //if(totalSegments<0) {cout<<"GOT TOTAL SEGMENT LESS THAN 0"<<endl; exit(0);}
    return totalSegments-1;
}

type_t BPlusTree::findCardinality(leaf *l, PMA *obj){
    type_t total = 0;
    for(int i=0; i<l->childCount; i++){
        /*
        if(l->childCount>Leaf_Degree){
            vector<node *> temp;
            temp.push_back(root);
            treeLevel = leafCount = 0;
            printTree(temp,1);
            exit(0);
        }
        */
        total += obj->cardinality[l->segNo[i]];
    }
    return total;
}

type_t BPlusTree::findCardinality(node *n, PMA *obj){
    type_t total = 0;
    for(int i=0; i < n->ptrCount; i++){
        /*
        if(n->ptrCount>Tree_Degree){
            vector<node *> temp;
            temp.push_back(root);
            treeLevel = leafCount = 0;
            printTree(temp, 1);
            exit(0);
        }
        */
        if(n->nodeLeaf) total += findCardinality((leaf *)n->child_ptr[i], obj);
        else total += findCardinality(n->child_ptr[i], obj);
    }
    return total;
}

void BPlusTree::printAllElements(PMA *obj){
    int totalElements = 0;
    leaf *leaf;
    for(leaf = leftmostLeaf(root); leaf != NULL; leaf = leaf->nextLeaf){
        for(int i = 0; i<leaf->childCount; i++){
            int segNo = leaf->segNo[i];
            type_t *key = obj->key_chunks[segNo];
            type_t pBase = 0;
            for(type_t block = 0; block<obj->blocksInSegment; block++){
                cout <<" Bitmap: "<<obj->bitmap[segNo][block]<<" ";
                u_short bitpos = 1;
                for(int j = 0; j<JacobsonIndexSize; j++){
                    if(obj->bitmap[segNo][block] & bitpos)
                        cout << *(key+pBase+j) << " ";
                    else cout <<"0 ";
                    bitpos = bitpos << 1;
                }
                pBase += JacobsonIndexSize;
                if(pBase>obj->lastElementPos[segNo]) break;
            }
            totalElements += obj->cardinality[segNo];
            cout<<"SIGMENT NO: "<<segNo<<" CARDINALITY: "<<obj->cardinality[segNo]<<" LAST Position: "<<obj->lastElementPos[segNo]<<endl;
        }
    }
    cout<<"Total element inserted in the PMA: "<<totalElements<<" Total Segments: "<<obj->totalSegments<<endl;
}

void BPlusTree::printTree(vector<Node *> nodes, int level){
    if(nodes[0] == NULL) return;
    if(nodes[0]->nodeLeaf){ //Create list of leaf child nodes from the current list
        cout<<"Printing level: "<<level<<endl;
        vector<Leaf *> list;
        for(u_int i=0; i<nodes.size(); i++){
            Node *temp = nodes[i];
            cout<<"Child:- "<<temp->ptrCount<<"::";
            for(int j=0; j<temp->ptrCount; j++){
                cout<<" "<<temp->child_ptr[j];
                if(j<temp->ptrCount-1) cout<<" "<<temp->key[j];
                list.push_back((Leaf *)temp->child_ptr[j]);
            }
            cout<<" || ";
        }
        cout<<endl;
        if(treeLevel == 0) treeLevel = level;
        else if(treeLevel != level){
            cout<<"Leaf Should have been at level: "<<treeLevel<<" Got at Level: "<<level<<endl;
        }
        leafCount++;
        printTree(list, level+1);
    }else{
        cout<<"Printing level: "<<level<<endl;
        vector<Node *> list;
        for(u_int i=0; i<nodes.size(); i++){
            Node *temp = nodes[i];
            cout<<"Child:- "<<temp->ptrCount<<"::";
            for(int j=0; j<temp->ptrCount; j++){
                cout<<" "<<temp->child_ptr[j];
                if(j<temp->ptrCount-1) cout<<" "<<temp->key[j];
                list.push_back(temp->child_ptr[j]);
            }
            cout<<" || ";
        }
        cout<<endl;
        printTree(list, level+1);
    }
}

void BPlusTree::printTree(vector<Leaf *> nodes, int level){
    cout<<"Printing level: "<<level<<endl;
    for(u_int i=0; i<nodes.size(); i++){
        Leaf *temp = nodes[i];
        cout<<"child:- "<<temp->childCount<<"::";
        for(int j=0; j<temp->childCount; j++){
            cout<<" "<<temp->segNo[j];
            if(j<temp->childCount-1) cout<<" "<<temp->key[j];
        }
        cout<<" || ";
    }
    cout<<endl;
}
