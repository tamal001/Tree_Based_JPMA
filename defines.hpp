#ifndef DEFINES_HPP_
#define DEFINES_HPP_

//#define SEGMENT_SIZE 2097152
#define CHUNK_SIZE 2097152UL
//#define CHUNK_SIZE 262144
//#define SEGMENT_SIZE 16384//1024//16384//32768

#define SEGMENT_SIZE 1024
//#define TOTAL_SEGMENTS 65536

#define type_t int64_t

#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)

#define PROTECTION (PROT_READ | PROT_WRITE)

#ifndef MAP_HUGETLB
#define MAP_HUGETLB 0x40000 /* arch specific */
#endif

//1 for mmap, 2 for malloc
#ifndef Allocation_type
#define Allocation_type 2
#endif

#define Tree_Degree 4
#define Leaf_Degree 5
#define MaxLevel 65

#ifdef __ia64__
#define ADDR (void *)(0x8000000000000000UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_FIXED)
#else
#define ADDR (void *)(0x0UL)
//#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB)
#define FLAGS (MAP_SHARED |MAP_ANONYMOUS | MAP_HUGETLB)
#endif

#define JacobsonIndexSize 16
#define JacobsonIndexCount 65536
//#define ProbLimit 32
#define ProbLimit (SEGMENT_SIZE * 10)/800 //Sizeof(type_t) = 8
#define MaxGap 3

#endif