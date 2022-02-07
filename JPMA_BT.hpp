#ifndef JPMA_HPP_
#define JPMA_HPP_

#include <vector>
#include <tuple>

#include "defines.hpp"
//#include "BPlusTree.hpp"
using namespace std;

class PMA;

class BPlusTree{
public:
    typedef struct Leaf{
        type_t key[Leaf_Degree-1];
        int segNo[Leaf_Degree];
        short childCount;
        Leaf *nextLeaf;
        Leaf() : childCount(0), nextLeaf(NULL) {}
    }leaf;

    //No node should have a combination of child of leaf and node
    typedef struct Node{
        type_t key[Tree_Degree-1];
        Node *child_ptr[Tree_Degree];
        bool nodeLeaf; //Last non-leaf node has value true
        short ptrCount;
        Node() : ptrCount(0), nodeLeaf(false){}
    }node;
    node *root;
    double level[100];
    //int maxElementInSegment;

    BPlusTree(PMA *obj);
    leaf* findLeaf(type_t search_key);
    int searchSegment(type_t search_key);
    void insertInTree(int chunkNo, type_t search_key, PMA *obj);
    void reinsertInTree(vector<int> &segments, type_t cardi, PMA *obj);
    void insert_in_parent(void *left, type_t search_key, void *right, type_t key_for_leaf);
    node * findParent(void *n, type_t key_parent);

    void calculateThreshold();
    void listSegments(vector<int> &segments, node *parent);
    type_t findCardinality(BPlusTree::leaf *l, PMA *obj);
    type_t findCardinality(BPlusTree::node *n, PMA *obj);
    void redistributeInsert(int segment, type_t Skey, PMA *obj);

    leaf* leftmostLeaf(node *root);
    leaf* rightmostLeaf(node *root);
    void deleteNode(node *parent);
    void deleteLeaf(leaf *l, type_t SKey);
    void printAllElements(PMA *obj);
    void printTree(vector<Node *> nodes, int level);
    void printTree(vector<Leaf *> nodes, int level);
};

class PMA{
public:
    type_t MaxThreshold;
    type_t MinThreshold;
    vector<type_t *> key_chunks;
    vector<type_t *> value_chunks;
    vector<type_t> smallest;
    vector<type_t> lastElementPos;
    vector<int> cardinality;
    int totalSegments;
    int elementsInSegment;
    u_char NonZeroEntries[JacobsonIndexCount][JacobsonIndexSize+1];
    vector<vector<u_short>> bitmap;
    BPlusTree *tree;
    type_t lastValidPos;             //Last accessible slot in each segment
    int freeSegmentCount;
    type_t blocksInSegment;
    vector<type_t *> freeKeySegmentBuffer;
    vector<type_t *> freeValueSegmentBuffer;
    int redisInsCount = 0, redisUpCount = 0;
    vector<type_t *> cleanSegments;

    PMA();
    ~PMA();

    //Library functions
    bool insert(type_t key, type_t value, int count= 0);
    bool remove(type_t key);
    bool lookup(type_t key);
    type_t range_sum(type_t startKey, type_t endKey);

    //Support functions
    int searchSegment(type_t key);
    tuple<type_t *, type_t *> getSegment();
    void preCalculateJacobson();
    void insertInPosition(type_t position, int targetSegment, type_t key, type_t value);
    bool backSearchInsert(type_t position, type_t key, type_t value, int targetSegment, int count);
    bool insertForward(type_t position, type_t key, type_t value, int targetSegment, int count); //Extra
    bool insertBackward(type_t position, type_t key, type_t value, int targetSegment, int count); //Extra
    bool insertAfterLast(type_t position, type_t key, type_t value, int targetSegment, type_t foundKey, int count);
    void deleteInPosition(type_t position, int targetSegment, type_t key);
    void deleteSegment(int targetSegment);
    type_t findLocation(type_t key, int targetSegment);
    type_t findLocation1(type_t key, int targetSegment);
    type_t findLocation2(type_t key, int targetSegment);
    //void redistribute(int targetSegment, redistribution type);
    //void redistributeWithPrev(int targetSegment);
    //void redistributeWithNext(int targetSegment);
    int redistributeWithDividing(int targetSegment);
    void swapElements(type_t targetSegment, type_t position, type_t adjust);

    //Testing functions
    void printStat();
    void printAllElements();
    void printSegElements(int targetSegment);
};

#endif