#ifndef GLOBAL_H
#define GLOBAL_H

// global_var_dec.h
// declaration of global variables
// Yun Dong 

extern char description_name[500];
extern float spacing;
extern float cw;
extern float cl;
extern float pitch;
extern float pitch_half;
extern float py_edge_low, py_edge_up, pz_edge_low, pz_edge_up;
extern float kvy_edge_low[100], kvy_edge_up[100], kvz_edge_low[100], kvz_edge_up[100];
extern float kernel_corner[100][2];
extern float kvx_edge_low, kvx_edge_up;
extern float vx, vy, vz;
extern float mu_lso;

extern int x_plane;
extern int layer;
extern int ncy;
extern int ncz;
extern int n_crsl2pxl;
extern int nkernel;
extern int nimgy;
extern int nimgz;
extern int nimgx;
extern int ncrow;
extern int nccol;
extern int ncy_k;
extern int ncz_k;
extern int ncrow_k;
extern int nccol_k;
extern int ngrid;

extern char kernel_name[100][10];
// maybe unnecessary
extern float eps;

#endif
