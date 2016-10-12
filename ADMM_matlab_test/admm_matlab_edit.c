#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#define _GNU_SOURCE
#include <getopt.h>
#define OUTFILE1 "DHAPET_original_forward.DAT"
#define OUTFILE2 "DHAPET_original_backward.DAT"
#include "global_var_dec.h"
#include "psf_analytic.h"

#define min(a,b)(((a)>(b))?(a):(b))
#define max(a,b)(((a)<(b))?(a):(b))

float s_time()
{
    return (float)clock()/CLOCKS_PER_SEC;
}

int nboundary(int cy, int cz)
{
    if( cy < ncy_k && cy >= 0 &&
            cz < ncz_k && cz >= 0    )
        return 1;
    else return 0;
}

int boundary(int cy, int cz)
{
    if( cy < ncy && cy >= 0 &&
            cz < ncz && cz >= 0    )
        return 1;
    else return 0;
}

int boundary_z(int cz)
{
    if(cz < ncz && cz >= 0) return 1;
    else return 0;
}

int boundary_y(int cy)
{
    if(cy < ncy && cy >= 0) return 1;
    else return 0;
}

void config()
{
    nkernel=16;
    nimgy = 4 * ncy;
    nimgz = 4 * ncz;
    ncrow = ncy * ncz;
    nccol = ncy * ncz;
    ncy_k = 2 * ncz;
    ncz_k = 2 * ncz;
    ncrow_k = ncy_k * ncz_k;
    nccol_k = ncrow_k;
    
    printf("  ncy=%d, ncz=%d\n",ncy,ncz);
    printf("  pitch=%.2f\n",pitch);
    printf("  nimgy=%d, nimgz=%d, nimgx=%d\n",nimgy,nimgz,nimgx);
    printf("  ngrid=%d\n",ngrid);
    printf("  spacing=%2.0fmm\n",spacing);
    printf("  index of xplane=%d\n",x_plane);
    //printf("  the object data will be saved with description: %s\n",description_name);
    py_edge_low = -ncy_k * cw / 2.;
    py_edge_up  =  ncy_k * cw / 2.;
    pz_edge_low = -ncz_k * cw / 2.;
    pz_edge_up  =  ncz_k * cw / 2.;
    vx = spacing / nimgx;
    vy = cw / (4);
    vz = cw / (4);
    
    kvx_edge_low = vx * (x_plane - 0.5);
    kvx_edge_up  = vx * (x_plane + 0.5);
    
}


//////////////////////////////////////////////////
//      Normalize 
//        Input : pNorm_LOR
//        Output: 
//////////////////////////////////////////////////
void norm_LOR(float* pNorm_LOR)
{
    int jlayer1=0;
    int jlayer2=0;
    
    int fx_plane = x_plane +60;
    TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));
    
    int i,j;
    
    FILE *inTAB,*inINCRE,*inSparse;
    
    
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            char name[500];
            //char kernel_name[100][10];
            sprintf(kernel_name[i*4+j],"y%dz%d",j,i);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
            inSparse = fopen(name,"rb");
            fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);
            
            int LOR_initial=pTable[nccol_k].bias;
            
            int iz,jy,kcz,kcy;
            for(iz=0;iz<ncz;iz++)
            {
                
                int shift_z= (iz*4-ncz*2)/4;
                for(jy=0;jy<ncy;jy++)
                {
                    
                    int shift_y = (jy*4-ncy*2)/4;
                    for(kcz=0;kcz<ncz;kcz++)
                    {
                        int cz1=kcz+52-shift_z;
                        for(kcy=0;kcy<ncy;kcy++)
                        {
                            int cy1=kcy+68-shift_y;
                            if(nboundary(cy1,cz1))
                            {
                                int cy2min=pTable[cz1*ncy_k+cy1].cymin+shift_y;
                                int cz2min=pTable[cz1*ncy_k+cy1].czmin+shift_z;
                                int c2total=pTable[cz1*ncy_k+cy1+1].bias;
                                //if(c2total>0) printf("c2total=%d\n",c2total);
                                // if(pTable[cz1*ncy_k/4+cy1].bias>0) printf("ktotal=%d\n",pTable[cz1*ncy_k/4+cy1].bias);
                                int ktotal;
                                for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
                                {//printf("123\n");
                                    int cy2=pIncre[ktotal]._y+cy2min;
                                    int cz2=pIncre[ktotal]._z+cz2min;
                                    //printf("cy2= %d cz2= %d ",cy2,cz2);
                                    
                                    if(boundary(cy2-68,cz2-52))
                                    {
                                        pNorm_LOR[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]+=pSparse[ktotal];
                                        // printf(" p: %f \n",phantom[59*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            fclose(inTAB);
            fclose(inINCRE);
            fclose(inSparse);
            free(pIncre);
            free(pSparse);
        }
    }
    free(pTable);
    
}

void norm_voxel(float* pNorm_voxel)
{
    int rx_plane = x_plane +60;
    int jlayer1=0;
    int jlayer2=0;
    TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));
    
    int i,j;
    
    FILE *inTAB,*inINCRE,*inSparse;
    
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            char name[500];
            //char kernel_name[100][10];
            sprintf(kernel_name[i*4+j],"y%dz%d",j,i);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
            inSparse = fopen(name,"rb");
            fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);
            
            int iz,jy,kcz,kcy;
            int LOR_initial=pTable[nccol_k].bias;
            for(iz=0;iz<ncz;iz++)
            {
                
                int shift_z = (iz*4-ncz*2)/4;
                for(jy=0;jy<ncy;jy++)
                {
                    
                    int shift_y= (jy*4-ncy*2)/4;
                    
                    for(kcz=0;kcz<ncz;kcz++)
                    {
                        int cz1=kcz+52-shift_z;
                        for(kcy=0;kcy<ncy;kcy++)
                        {
                            int cy1=kcy+68-shift_y;
                            if(nboundary(cy1,cz1))
                            {
                                int cy2min=pTable[cz1*ncy_k+cy1].cymin+shift_y;
                                int cz2min=pTable[cz1*ncy_k+cy1].czmin+shift_z;
                                int c2total=pTable[cz1*ncy_k+cy1+1].bias;
                                int ktotal;
                                for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
                                {
                                    int cy2=pIncre[ktotal]._y+cy2min;
                                    int cz2=pIncre[ktotal]._z+cz2min;
                                    if(boundary(cy2-68,cz2-52))
                                    {
                                        
                                        pNorm_voxel[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pSparse[ktotal];
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
            fclose(inTAB);
            fclose(inINCRE);
            fclose(inSparse);
            free(pIncre);
            free(pSparse);
        }
    }
    free(pTable);
}


void forward(float* pFwd,float* phantom,float* pNorm_LOR)
{
    int jlayer1=0;
    int jlayer2=0;
    
    int fx_plane = x_plane +60;
    TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));
    
    int i,j;
    
    FILE *inTAB,*inINCRE,*inSparse;
    
    
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            char name[500];
            //char kernel_name[100][10];
            sprintf(kernel_name[i*4+j],"y%dz%d",j,i);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
            inSparse = fopen(name,"rb");
            fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);
            
            int LOR_initial=pTable[nccol_k].bias;
            
            int iz,jy,kcz,kcy;
            for(iz=0;iz<ncz;iz++)
            {
                
                int shift_z= (iz*4-ncz*2)/4;
                for(jy=0;jy<ncy;jy++)
                {
                    
                    int shift_y = (jy*4-ncy*2)/4;
                    for(kcz=0;kcz<ncz;kcz++)
                    {
                        int cz1=kcz+52-shift_z;
                        for(kcy=0;kcy<ncy;kcy++)
                        {
                            int cy1=kcy+68-shift_y;
                            if(nboundary(cy1,cz1))
                            {
                                int cy2min=pTable[cz1*ncy_k+cy1].cymin+shift_y;
                                int cz2min=pTable[cz1*ncy_k+cy1].czmin+shift_z;
                                int c2total=pTable[cz1*ncy_k+cy1+1].bias;
                                //if(c2total>0) printf("c2total=%d\n",c2total);
                                // if(pTable[cz1*ncy_k/4+cy1].bias>0) printf("ktotal=%d\n",pTable[cz1*ncy_k/4+cy1].bias);
                                int ktotal;
                                for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
                                {//printf("123\n");
                                    int cy2=pIncre[ktotal]._y+cy2min;
                                    int cz2=pIncre[ktotal]._z+cz2min;
                                    //printf("cy2= %d cz2= %d ",cy2,cz2);
                                    
                                    if(boundary(cy2-68,cz2-52))
                                    {
                                        pFwd[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]+=pSparse[ktotal]*phantom[fx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]/pNorm_LOR[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)];///(float)LOR_initial;
                                        // printf(" p: %f \n",phantom[59*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            fclose(inTAB);
            fclose(inINCRE);
            fclose(inSparse);
            free(pIncre);
            free(pSparse);
        }
    }
    free(pTable);
    
}

void ratio(float *pFW, float *pHist, float *pRatio)
{
    int i;
    for(i = 0; i < ncrow * nccol; i++)
    {
        if(pHist[i] >= 0)
        {
            pRatio[i] =  abs(pFW[i]-pHist[i]);
            
            pRatio[i]=2*pRatio[i];
        }
        
    }
    
}
/////////////Edit by PeiHsiu/////////////////
// void norm(float *sol, float* norm2,int ite)
// {
//     int i;
    
//     for(i=0;i<nimgx*nimgy*nimgz;i++)
//     {
//         norm2[ite] = sqrt( norm2[ite]*norm2[ite] + sol[i]*sol[i] );
//     }
    
// }


//////////////////////////////////////////////////Edit by PeiHsiu
//      2-Dimentional TV norm
//        Input : sol
//        Output: normtv
//////////////////////////////////////////////////
// void norm_tv(float *sol, float *normtv, float* ite)
// {
//     //////////////////start gradient_op/////////////////////////////////
//     int temp[nimgy*nimgz]; int dx[nimgy*nimgz]; int dy[nimgy*nimgz]; int y;
//     int x_n[nimgy*nimgz+1];
//     int i;  int k;
//     int rx_plane = x_plane + 60;
//     for (i=0;i<nimgy*nimgz;i++)
//     {
//         if (i%nimgy==nimgy-1)
//             dx[i] = 0;
//         else
//             dx[i] = x_n[rx_plane*nimgy*nimgz+i+1] - x_n[rx_plane*nimgy*nimgz+i];
//     }
    
//     for (i=0;i<nimgz;i++){
//         for (k=0;k<nimgy;k++){
//             if (i==nimgz-1)
//                 dy[i*nimgy+k] = 0;
//             else
//                 dy[i*nimgy+k] = x_n[rx_plane*nimgy*nimgz+(i+1)*nimgy+k] - x_n[rx_plane*nimgy*nimgz+i*nimgy+k];            
//         }
//     }

//     for (i=0;i<nimgy*nimgz;i++){
//         temp[i] = sqrt( dx[i]*dx[i] + dy[i]*dy[i] );
//         y = sum(sum(temp,1),2);
//     }
    
//     // int i;
//     for(i=0;i<nimgx-1;i++){
//         for(j=0;j<nimgy-1;j++){
//             for(k=0;k<nimgz;k++){
//                 dx[i][j][k] = x[i+1][j][k] - x[i][j][k] ;//, zeros(1,sizeof(x,2),sizeof(x,3))];
//                 dy[i][j][k] = x[i][j+1][k] - x[i][j][k] ;//, zeros(sizeof(x,1), 1,sizeof(x,3))];
//             }
//         }
//     }

//     temp[ite] = sqrt( dx[ite]*dx[ite] + dy[ite]*dy[ite] );
//     y = sum(sum(temp,1),2);
//     }
    
// }

/////////////End of Edit///////////////////////

void cost_function(float * pFW, float *pHist, float* norm2, float *costfunction, int ite, float gamma)
{
    int i;
    for(i = 0; i < ncrow * nccol; i++)
    {
        
        costfunction[ite] =  costfunction[ite] + pow((pFW[i]-pHist[i]),2);
        
        
    }
    printf("first term[%d] = %f\n", ite, costfunction[ite]);
    printf("second term = %f\n",gamma*sqrt(norm2[ite]));
    costfunction[ite] =  costfunction[ite] + gamma *sqrt(norm2[ite]);
    printf("cost function[%d] = %f\n", ite, costfunction[ite]);
}

void backward(float* pRatio, float* pRecon, float* pNorm_voxel)
{
    int rx_plane = x_plane +60;
    int jlayer1=0;
    int jlayer2=0;
    TABLE_ZIP* pTable = (TABLE_ZIP*) calloc(nccol_k+1, sizeof(TABLE_ZIP));
    
    int i,j;
    
    FILE *inTAB,*inINCRE,*inSparse;
    
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            char name[500];
            //char kernel_name[100][10];
            sprintf(kernel_name[i*4+j],"y%dz%d",j,i);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../hsuyueh/two_level_EM/EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
            inSparse = fopen(name,"rb");
            fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);
            
            int iz,jy,kcz,kcy;
            int LOR_initial=pTable[nccol_k].bias;
            for(iz=0;iz<ncz;iz++)
            {
                
                int shift_z = (iz*4-ncz*2)/4;
                for(jy=0;jy<ncy;jy++)
                {
                    
                    int shift_y= (jy*4-ncy*2)/4;
                    
                    for(kcz=0;kcz<ncz;kcz++)
                    {
                        int cz1=kcz+52-shift_z;
                        for(kcy=0;kcy<ncy;kcy++)
                        {
                            int cy1=kcy+68-shift_y;
                            if(nboundary(cy1,cz1))
                            {
                                int cy2min=pTable[cz1*ncy_k+cy1].cymin+shift_y;
                                int cz2min=pTable[cz1*ncy_k+cy1].czmin+shift_z;
                                int c2total=pTable[cz1*ncy_k+cy1+1].bias;
                                int ktotal;
                                for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
                                {
                                    int cy2=pIncre[ktotal]._y+cy2min;
                                    int cz2=pIncre[ktotal]._z+cz2min;
                                    if(boundary(cy2-68,cz2-52))
                                    {
                                        pRecon[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pRatio[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]*pSparse[ktotal]/pNorm_voxel[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j];
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
            fclose(inTAB);
            fclose(inINCRE);
            fclose(inSparse);
            free(pIncre);
            free(pSparse);
        }
    }
    free(pTable);
}



void proxtv(float* x_n,float* u_n, float* r, float *s, float* r2, float* s2,float* dx, float* dy,float* pold, float* qold, float* p, float* q, float* temp_s, float* temp_r, float* sum_rs,float gamma, float* Y_n)

{
    float weights;
    
    float weight_x=1,weight_y=1,mt=1,t,told=1;
    int rx_plane = x_plane +60;
    int i,k;
 /////////////////////proxTV iteration//////////////////////////////
    
    //printf("%d",iteration);
    for (i=0;i<nimgy*nimgz;i++)
    {
        if (i%nimgy==0)
        {	s2[i]=s[i];}
        else if (i%nimgy==nimgy-1)
        {	s2[i]=-s[i-1];}
        else
        {s2[i]=s[i]-s[i-1];}
        
        s2[i]=s2[i]*weight_x;
        
    }
    
    for (i=0;i<nimgz;i++)
    {
        for (k=0;k<nimgy;k++)
        {
            if (i==0)
            {r2[i*nimgy+k]=r[i*nimgy+k];}
            else if (i==nimgz-1)
            {r2[i*nimgy+k]=-1*r[(i-1)*(nimgy)+k];}
            else
            {r2[i*nimgy+k]=r[i*nimgy+k]-r[(i-1)*nimgy+k];}
            
            r2[i]=r2[i]*weight_y;
            
        }
    }
     /////////////current solution//////////////
    for (i=0;i<nimgy*nimgz;i++)
    {
        x_n[rx_plane*nimgy*nimgz+i]=Y_n[i]-gamma*(r2[i]+s2[i]);
        //sol[rx_plane*nimgy*nimgz+i]=Y_n[i]-gamma*(r2[i]+s2[i]);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////
    //update divergence vectors and project
    //////////////////start gradient_op/////////////////////////////////
    for (i=0;i<nimgy*nimgz;i++)
    {
        if (i%nimgy==nimgy-1)
            dx[i]=0;
        else
            dx[i]=x_n[rx_plane*nimgy*nimgz+i+1]-x_n[rx_plane*nimgy*nimgz+i];
    }
    
    for (i=0;i<nimgz;i++)
    {
        for (k=0;k<nimgy;k++)
        {
            if (i==nimgz-1)
                dy[i*nimgy+k]=0;
            else
                dy[i*nimgy+k]=x_n[rx_plane*nimgy*nimgz+(i+1)*nimgy+k]-x_n[rx_plane*nimgy*nimgz+i*nimgy+k];
            
        }
    }
 ////////////////////////////////////////////////////////////////////////////////////////////
    float max=0;
    float weight;
    for (i=0;i<nimgy*nimgz;i++)
    {
        r[i]=r[i]-(1/(8*gamma)/(mt*mt)*dy[i]);
        s[i]=s[i]-(1/(8*gamma)/(mt*mt)*dx[i]);
        temp_r[i]=abs(r[i]);
        temp_s[i]=abs(s[i]);
        sum_rs[i]=temp_r[i]*temp_r[i]+temp_s[i]*temp_s[i];
        sum_rs[i]=sqrt(sum_rs[i]);
        
        if (sum_rs[i]>=max)
        {
            max=sum_rs[i];
        }

        weights=max;

        if (max<1)
            weights=1;
    }
    
    for (i=0;i<nimgy*nimgz;i++)
    {
        p[i]=r[i]/weights;
        q[i]=s[i]/weights;  
    }
    
    // FISTA update
    t=(1+sqrt(4*told*told))/2;
    
    for (i=0;i<nimgy*nimgz;i++)
    {
        r[i]=p[i]+((told-1)/t)*(p[i]-pold[i]);
        pold[i]=p[i];
        s[i]=q[i]+((told-1)/t)*(q[i]-qold[i]);
        qold[i]=q[i];
    }
    
    told=t;
}



/////////////////////proj_B2 iteration//////////////////////////////by PeiHsiu
void proj_b2(float* s2, float* s, float* r2, float* r, float* f2, float* u, float* x_n, float* res, float* temp, float* x, float* Y_n, float* sctemp, float* sol)
{
    float told = 1;
    int i, j, k;
    //float res, u = 0,
    float v, ratio, t;
    //float fminf(float temp);
    float nu = 1;

    float weight_x=1, weight_y=1, mt=1;
    float epsilon[nimgy*nimgz];
    int rx_plane = x_plane + 60;
    float* temp_min = (float*)calloc(nimgy*nimgz,sizeof(float));
    /////From Prox_TV//////

    for (i = 0; i < nimgy*nimgz; i++){
        if (i%nimgy==0)
        {   s2[i] = s[i];}
        else if (i%nimgy==nimgy-1)
        {   s2[i] = -s[i-1];}
        else
        {   s2[i] = s[i]-s[i-1];}
        
        s2[i] = s2[i]*weight_x;
    }
    
    for (i = 0; i < nimgz; i++){
        for (k = 0; k < nimgy; k++){
            if (i==0)
            {r2[i*nimgy+k] = r[i*nimgy+k];}
            else if (i==nimgz-1)
            {r2[i*nimgy+k] = -1*r[(i-1)*(nimgy)+k];}
            else
            {r2[i*nimgy+k] = r[i*nimgy+k]-r[(i-1)*nimgy+k];}
            
            r2[i] = r2[i]*weight_y;
        }
    }

    ///////////Projection//////////////
    //Tight Frame Case
    float norm2[nimgy*nimgz];
    for (i=0;i<nimgy*nimgz;i++)
    {   
        temp[i] = x_n[i] - Y_n[i];
        for (i=0;i<nimgy*nimgz;i++)
        {
            temp_min = epsilon[i]/norm(temp[i],norm2[i]);
        }
        sctemp[i] = temp[i] * min(temp_min);   //Scaling
        sol[i] = x[i] + 1/nu * (sctemp[i] - temp[i]);
        temp[i] = sol[i];
    }




    /*
    //Non Tight Frame Case
    // Projection onto the L2-ball
    
    // Update number of iteration
    //iter = iter + 1;

    for (i = 0; i < nimgy*nimgz; i++)
    {
        // Residual
        res[i] = x[i] - y[i];
        norm_res = norm(res,2);
        
        // Scaling for the projection
        res[i] = u[i] * nu + res[i];
        norm_proj = norm(res,2);
    }
        
    ratio = min(1, epsilon/norm_proj);
        
    t = (1 + sqrt(1 + 4*told^2))/2;		//FISTA timestop
    for (i = 0; i < nimgy*nimgz; i++)
    {
        u[i] = v[i];
        v[i] = 1/nu * (res[i] - res[i]*ratio);
        u[i] = v[i] + (told - 1)/t * (v[i] - u[i]);
    }

    // Update timestop
    told = t;
        
    // Current estimate
    sol = x - At(u);
    }
    */ 
}




int main()
{
    config();
    //pHist_1:nonTOF
    float* pHist_1= (float*) calloc(nccol*ncrow,sizeof(float));
    //pHist_2:TOF
    float* pHist_2= (float*) calloc(nccol*ncrow,sizeof(float));
    
    float* pFwd_1= (float*) calloc(nccol*ncrow,sizeof(float));
    float* pRatio_1= (float*) calloc(nccol*ncrow,sizeof(float));
    float* phantom_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pBwd_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* sol_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* x_n_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* u_n_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pRecon_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    
    float* pFwd_2= (float*) calloc(nccol*ncrow,sizeof(float));
    float* pRatio_2= (float*) calloc(nccol*ncrow,sizeof(float));
    float* phantom_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pBwd_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* sol_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* x_n_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* u_n_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pRecon_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    
    float* pNorm_voxel_1= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pNorm_LOR_1= (float*) calloc(nccol*ncrow,sizeof(float));
    float* costfunction_1= (float*) calloc(100,sizeof(float));
    float* norm2_1=(float*) calloc(100,sizeof(float));
    
    float* pNorm_voxel_2= (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pNorm_LOR_2= (float*) calloc(nccol*ncrow,sizeof(float));
    float* costfunction_2= (float*) calloc(100,sizeof(float));
    float* norm2_2=(float*) calloc(100,sizeof(float));
    
    float * Y_n_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * r_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * s_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * dx_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * dy_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * p_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * q_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * pold_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * qold_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_r_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_s_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * sum_rs_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * r2_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * s2_1=(float*)calloc(nimgy*nimgz,sizeof(float));
    
    float * Y_n_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * r_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * s_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * dx_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * dy_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * p_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * q_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * pold_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * qold_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_r_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_s_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * sum_rs_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * r2_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    float * s2_2=(float*)calloc(nimgy*nimgz,sizeof(float));
    
    ////////////////////////////Edit by PeiHsiu/////////////////////////
    float * s1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * f2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * u = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * x = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * res = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp = (float*)calloc(nimgy*nimgz,sizeof(float));
    // float * dx_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    // float * dy_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * sctemp = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * sol = (float*)calloc(nimgy*nimgz,sizeof(float));
    ////////////////////////////////////////////////////////////////////


    float tn=1,tn1;
    int i,j,k;
    float gamma=0.01;
    int iteration=2;
    int rx_plane = x_plane +60;
    int a1=-675; int a2=-4.95;
    
    int ix,iy,iz,kk,ss;
    
    memset(pFwd_1,0, ncrow*nccol * sizeof(float));
    memset(costfunction_1,0, 100*sizeof(float));
    memset(norm2_1,0, 100*sizeof(float));
    memset(pFwd_2,0, ncrow*nccol * sizeof(float));
    memset(costfunction_2,0, 100*sizeof(float));
    memset(norm2_2,0, 100*sizeof(float));
    //printf("stop 01\n");
    
    FILE *iHist_1;
    iHist_1=fopen("../histogram/histogram_nonTOF_new.dat","rb");
    fread(pHist_1,sizeof(float),nccol*ncrow,iHist_1);
    fclose(iHist_1);
    
    FILE *iHist_2;
    iHist_2=fopen("../histogram/histogram_TOF_new.dat","rb");
    fread(pHist_2,sizeof(float),nccol*ncrow,iHist_2);
    fclose(iHist_2);
    // printf("stop 02\n");
    
    
    for(ix=0;ix<120;ix++)
    {
        for(i=0;i<nimgy*nimgz;i++)
        {
            pBwd_1[ix*nimgz*nimgy+i] = 0.0;
            x_n_1[ix*nimgz*nimgy+i] = pHist_1[ix*nimgz*nimgy+i];
            sol_1[ix*nimgz*nimgy+i] = 0.0;
            u_n_1[ix*nimgz*nimgy+i] = 0.0;

            pBwd_2[ix*nimgz*nimgy+i] = 0.0;
            x_n_2[ix*nimgz*nimgy+i] = pHist_2[ix*nimgz*nimgy+i];
            sol_2[ix*nimgz*nimgy+i] = 0.0;
            u_n_2[ix*nimgz*nimgy+i] = 0.0;
        }
    }
    
    
    
 /////////////////////////////////start ADMM////////////////////////////////////////////
    int ite;
    float T[2];
    char out_name[300];
    
    T[0]=s_time();
    memset(pFwd_1,0,nccol*ncrow*sizeof(float));
    memset(pRecon_1,0,nimgx*nimgy*nimgz*sizeof(float));
    memset(pNorm_voxel_1,0,nimgx*nimgy*nimgz*sizeof(float));
    memset(pNorm_LOR_1,0,nccol*ncrow*sizeof(float));
    
    memset(pFwd_2,0,nccol*ncrow*sizeof(float));
    memset(pRecon_2,0,nimgx*nimgy*nimgz*sizeof(float));
    memset(pNorm_voxel_2,0,nimgx*nimgy*nimgz*sizeof(float));
    memset(pNorm_LOR_2,0,nccol*ncrow*sizeof(float));
    
    for(x_plane=0;x_plane<3;x_plane++)
    {
        x_plane=x_plane*(-20);
        norm_LOR(pNorm_LOR_1);
        norm_voxel(pNorm_voxel_1);
        norm_LOR(pNorm_LOR_2);
        norm_voxel(pNorm_voxel_2);
        x_plane=x_plane/(-20);
    }
    
    for(x_plane=0;x_plane<3;x_plane++)
    {
        x_plane=x_plane*(-20);
        //printf("start forward: %d\n",x_plane);
        forward(pFwd_1,pBwd_1,pNorm_LOR_1);
        forward(pFwd_2,pBwd_2,pNorm_LOR_2);
        x_plane=x_plane/(-20);
    }
    
    //printf("start ratio\n");
    ratio(pFwd_1,pHist_1,pRatio_1);
    ratio(pFwd_2,pHist_2,pRatio_2);
    
    for(x_plane=0;x_plane<3;x_plane++)
    {
        x_plane=x_plane*(-20);
        //printf("start backward: %d\n",x_plane);
        backward(pRatio_1,pRecon_1,pNorm_voxel_1);
        backward(pRatio_2,pRecon_2,pNorm_voxel_2);
        x_plane=x_plane/(-20);
    }
    

    /////////Edit by PeiHsiu//////////
    float curr_norm[nimgy][nimgz];
    int tau = 50;
    float x_0[nimgy];
    float norm2[nimgy*nimgz];
    float prev_norm[nimgy][nimgz];
    float y_old[nimgy*nimgz];
    float y_n[nimgy*nimgz];
    float z_n[nimgy*nimgz];

    for(i=0;i<nimgy*nimgz;i++){
        curr_norm[i] = tau * norm_tv(x_0[i]) + norm(x_0[i],norm2,i);
        prev_norm[i] = convergence_test(curr_norm[i]);
        y_old[i] = y_n[i];
        z_n[i] = zeros(sizeof(y_n));
    }

    float tol = 10e-4;
    float reldual = 1;

    while(reldual < tol)
    {
    for (ite=0;ite<iteration;ite++)
    {
        
     /////////////////u_n-gamma*f2.grad(u_n)///////////////////////////////
        for(x_plane=0;x_plane<3;x_plane++)
        {
            x_plane=x_plane*(-20);
            for (i=0;i<nimgy*nimgz;i++)
            {
                rx_plane = x_plane +60;
                Y_n_1[i]=pBwd_1[rx_plane*nimgy*nimgz+i]-gamma*pRecon_1[rx_plane*nimgy*nimgz+i];
                
                /////////////////initial/////////////////////
                r_1[i]=0;
                s_1[i]=0;
                pold_1[i]=r_1[i];
                qold_1[i]=s_1[i];
            }

            //Main loop

            //printf("iteration=%d\n",ite+1);
            proxtv(x_n_1,u_n_1,r_1,s_1,r2_1,s2_1,dx_1,dy_1,pold_1,qold_1,p_1,q_1,temp_s_1,temp_r_1,sum_rs_1,gamma,Y_n_1);
            
            proj_b2(s1, f2, r2, r, u, x, res, temp, dx_1, dy_1, sctemp, sol);

            float s_n[nimgy*nimgz];
            float x_n[nimgy*nimgz];

            for(i=0;i<nimgy*nimgz;i++){
                reldual[i] = norm(y_old[i] - y_n[i])/norm(y_n[i]);

                z_n[i] = z_n[i] + s_n[i] - y_n[i];  //Updates
                sol[i] = x_n[i];
                y_old[i] = y_n[i];
            }

            for(i=0;i<nimgy*nimgz;i++)
            {
                sol_1[rx_plane*nimgy*nimgz+i]=a1*x_n_1[rx_plane*nimgy*nimgz+i];
            }
            
            if(((ite+1)%5==0) || (ite<5))
            {
                if(x_plane==-40)
                {
                    forward(pFwd_1,sol_1,pNorm_LOR_1);
                    norm(sol_1,norm2_1,ite);
                    //printf("norm2=%f\n",norm2[ite]);
                    cost_function(pFwd_1,pHist_1,norm2_1,costfunction_1,ite,gamma);
                    
                    sprintf(out_name,"GPU_DHAPET_backward_ite%d_activity_r001_new.DAT",ite+1);
                    FILE *fp1;
                    if((fp1=fopen(out_name,"wb"))==NULL)
                    {
                        printf("Can not open file\n");
                        exit(0);
                    }
                    fwrite(sol_1,sizeof(float)*nimgx*nimgy*nimgz,1,fp1);
                    fclose(fp1);
                }
            }
            x_plane=x_plane/(-20);
        }
    }
    }   //End of while loop
    //Global stopping criterion
    for(i=0;i<nimgy*nimgz;i++){
        curr_norm[i] = tau * norm_tv(sol[i]) + norm(sol[i],norm2,ite);
    }










 //////////////////////////////////Attenuation/////////////////////////////////


    for (ite=0;ite<iteration;ite++)
    {
        
/////////////////u_n-gamma*f2.grad(u_n)///////////////////////////////
        for(x_plane=0;x_plane<3;x_plane++)
        {
            x_plane=x_plane*(-20);
            for (i=0;i<nimgy*nimgz;i++)
            {
                rx_plane = x_plane +60;
                Y_n_2[i]=pBwd_2[rx_plane*nimgy*nimgz+i]-gamma*pRecon_2[rx_plane*nimgy*nimgz+i];
/////////////////initial////////////////////////////////////////////
                r_2[i]=0;
                s_2[i]=0;
                pold_2[i]=r_2[i];
                qold_2[i]=s_2[i];
            }
            
            //printf("iteration=%d\n",ite+1);
            proxtv(x_n_2,u_n_2,r_2,s_2,r2_2,s2_2,dx_2,dy_2,pold_2,qold_2,p_2,q_2,temp_s_2,temp_r_2,sum_rs_2,gamma,Y_n_2);
            
            proj_b2(s1, f2, r2, r, u, x, res, temp, x, y, sctemp, sol);

            for(i=0;i<nimgy*nimgz;i++)
            {
                sol_2[rx_plane*nimgy*nimgz+i]=a2*x_n_2[rx_plane*nimgy*nimgz+i];
            }
            
            if(((ite+1)%5==0) || (ite<5))
            {
                if(x_plane==-40)
                {
                    forward(pFwd_2,sol_2,pNorm_LOR_2);
                    norm(sol_2,norm2_2,ite);
                    //printf("norm2_2=%f\n",norm2_2[ite]);
                    cost_function(pFwd_2,pHist_2,norm2_2,costfunction_2,ite,gamma);
                    
                    sprintf(out_name,"GPU_DHAPET_backward_ite%d_attenuation_r001_new.DAT",ite+1);
                    FILE *fp2;
                    if((fp2=fopen(out_name,"wb"))==NULL)
                    {
                        printf("Can not open file\n");
                        exit(0);
                    }
                    fwrite(sol_2,sizeof(float)*nimgx*nimgy*nimgz,1,fp2);
                    fclose(fp2);
                }
            }
            x_plane=x_plane/(-20);
        }
    }
    T[1]=s_time()-T[0];
    printf("time: %g sec \n",T[1]);
    
    
/////////////////////////////////stop em/////////////////////////////////////////////
    
    
    free(norm2_1);
    free(x_n_1);
    free(u_n_1);
    free(sol_1);
    free(phantom_1);
    free(pFwd_1);
    free(pBwd_1);
    free(pHist_1);
    free(pRatio_1);
    free(pRecon_1);
    free(pNorm_LOR_1);
    free(pNorm_voxel_1);
    free(costfunction_1);
    
    free(norm2_2);
    free(x_n_2);
    free(u_n_2);
    free(sol_2);
    free(phantom_2);
    free(pFwd_2);
    free(pBwd_2);
    free(pHist_2);
    free(pRatio_2);
    free(pRecon_2);
    free(pNorm_LOR_2);
    free(pNorm_voxel_2);
    free(costfunction_2);
    
    free(r_1);
    free(s_1);
    free(dx_1);
    free(dy_1);
    free(p_1);
    free(q_1);
    free(qold_1);
    free(pold_1);
    free(temp_r_1);
    free(temp_s_1);
    free(sum_rs_1);
    free(r2_1);
    free(s2_1);
    free(Y_n_1);
    
    free(r_2);
    free(s_2);
    free(dx_2);
    free(dy_2);
    free(p_2);
    free(q_2);
    free(qold_2);
    free(pold_2);
    free(temp_r_2);
    free(temp_s_2);
    free(sum_rs_2);
    free(r2_2);
    free(s2_2);
    free(Y_n_2);
}
