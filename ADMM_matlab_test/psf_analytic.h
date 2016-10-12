#ifndef PSF_H
#define PSF_H

// psf_analytic.h
// histogram save and read as phist[104(highest_dimension)][72][104][72(lowest_dimension)], or 104*72*104*72
// Yun Dong 

typedef struct
{
    int cyk1;
    int czk1;
    int cyk2;
    int czk2;
    int s;      // simulated voxel; vy=0, vz=0 (s=0); vy=1, vz=0(s=1); vy=1, vz=1(s=2)
} kernel;

typedef struct
{
    int cymin1;
    int cymax1;
    int czmin1;
    int czmax1;
    int bias1;
} table1;

typedef struct
{
    int cymin2;
    int cymax2;
    int czmin2;
    int czmax2;
    int bias2;
} table2;

typedef struct
{
    int cymin;
    int cymax;
    int czmin;
    int czmax;
    int bias;
} table;

typedef struct
{
    int cymin1;
    int czmin1;
    int bias1;
} TABLE_ZIP1;

typedef struct
{
    int cymin2;
    int czmin2;
    int bias2;
} TABLE_ZIP2;

typedef struct
{
    int cymin;
    int czmin;
    int bias;
} TABLE_ZIP;

typedef struct
{
    unsigned char _y;
    unsigned char _z;
} INCREMENT_ZIP1;

typedef struct
{
    unsigned char _y;
    unsigned char _z;
} INCREMENT_ZIP2;

typedef struct
{
    unsigned char _y;
    unsigned char _z;
} INCREMENT_ZIP;

typedef struct
{
    TABLE_ZIP1*** tab1;
    TABLE_ZIP2*** tab2;
    TABLE_ZIP*** tab;
    INCREMENT_ZIP1*** inc1;
    INCREMENT_ZIP2*** inc2;
    INCREMENT_ZIP*** inc;
    float*** spm;
    char** tag_mem_tab;
    char** tag_mem_inc;
    char** tag_mem_spm;
    int call_tag;
}MEMCALLOC;

typedef struct
{
    float alpha;
    char io;
}TRACING_CRY;

typedef struct
{
    int cy;
    int cz;
    float resp;
}RESPONSE_CRY;

typedef struct
{
    unsigned char layer1;
    unsigned char layer2;
    unsigned char cy1;
    unsigned char cz1;
    unsigned char cy2;
    unsigned char cz2;
    float co_resp;
}CO_RESPONSE;

// low_precision.c
float single_ray_tracing_kvoxel(float, float, float, float, int);
void tracing_kvoxel();

// high_precision.c

#endif
