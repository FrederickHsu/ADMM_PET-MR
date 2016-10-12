#include "global_var_dec.h"

// global_var_def.c: definition of global variables
// Yun Dong 

char description_name[500];

float spacing = 60;	//distance between two plate
float cw = 2.42;	//crystal dimention (y z direction)
float cl = 10;		//half of crystal length (x direction)
//float pitch = 0;
float pitch = 0.3;	//distance between two crystal
float pitch_half;
float py_edge_low, py_edge_up, pz_edge_low, pz_edge_up;
float kvy_edge_low[100], kvy_edge_up[100], kvz_edge_low[100], kvz_edge_up[100];
float kvx_edge_low, kvx_edge_up;
float kernel_corner[100][2];
float vx, vy, vz;
float mu_lso = 0.088; // LSO linear atenuation coefficients at 5llkeV , /mm

int x_plane = 0;
int layer = 2;
int ncy = 72 ; 
int ncz = 104; 
int n_crsl2pxl = 4;
int nkernel;	//number of simulating grid of object space

int ncrow;
int nccol;

int ncy_k;
int ncz_k;
int ncrow_k;
int nccol_k;

int ngrid = 4;

int nimgx = 120;
int nimgy;
int nimgz;
int imgybias; 
int imgzbias;

int cybias_k;
int czbias_k;

int vybias_k;
int vzbias_k;

float eps = 1.0e-6;

char kernel_name[100][10];

