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
    nkernel = 16;
    nimgy = 4 * ncy;
    nimgz = 4 * ncz;
    ncrow = ncy * ncz;
    nccol = ncy * ncz;
    ncy_k = 2 * ncz;
    ncz_k = 2 * ncz;
    ncrow_k = ncy_k * ncz_k;
    nccol_k = ncrow_k;
    
    printf("  ncy = %d, ncz = %d\n",ncy,ncz);
    printf("  pitch = %.2f\n",pitch);
    printf("  nimgy = %d, nimgz = %d, nimgx = %d\n",nimgy,nimgz,nimgx);
    printf("  ngrid = %d\n",ngrid);
    printf("  spacing = %2.0fmm\n",spacing);
    printf("  index of xplane = %d\n",x_plane);
    printf("  the object data will be saved with description: %s\n",description_name);
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

// Normalize LOR
void norm_LOR(float *pNorm_LOR)
{
    int jlayer1 = 0;
    int jlayer2 = 0;
    
    int fx_plane = x_plane + 60;
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
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
            inSparse = fopen(name,"rb");
            fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);
            
            int LOR_initial = pTable[nccol_k].bias;
            
            int iz,jy,kcz,kcy;
            for(iz=0;iz<ncz;iz++)
            {
                
                int shift_z = (iz*4-ncz*2)/4;
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
                                int cy2min = pTable[cz1*ncy_k+cy1].cymin + shift_y;
                                int cz2min = pTable[cz1*ncy_k+cy1].czmin + shift_z;
                                int c2total = pTable[cz1*ncy_k+cy1+1].bias;
                                //if(c2total>0) printf("c2total=%d\n",c2total);
                                // if(pTable[cz1*ncy_k/4+cy1].bias>0) printf("ktotal=%d\n",pTable[cz1*ncy_k/4+cy1].bias);
                                int ktotal;
                                for(ktotal=pTable[cz1*ncy_k+cy1].bias;ktotal<c2total;ktotal++)
                                {//printf("123\n");
                                    int cy2 = pIncre[ktotal]._y + cy2min;
                                    int cz2 = pIncre[ktotal]._z + cz2min;
                                    // printf("cy2= %d cz2= %d ",cy2,cz2);
                                    
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

// Normalize Voxel
void norm_voxel(float *pNorm_voxel)
{
    int rx_plane = x_plane +60;
    int jlayer1 = 0;
    int jlayer2 = 0;
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
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse= (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
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


void forward(float *pFwd, float *phantom, float *pNorm_LOR)
{
    int jlayer1 = 0;
    int jlayer2 = 0;
    
    int fx_plane = x_plane + 60;
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
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse = (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
            inSparse = fopen(name,"rb");
            fread(pSparse,sizeof(float),pTable[nccol_k].bias,inSparse);
            
            int LOR_initial = pTable[nccol_k].bias;
            
            int iz, jy, kcz, kcy;
            for(iz=0;iz<ncz;iz++)
            {
                
                int shift_z = (iz*4-ncz*2)/4;
                for(jy=0;jy<ncy;jy++)
                {
                    
                    int shift_y = (jy*4-ncy*2)/4;
                    for(kcz=0; kcz<ncz; kcz++)
                    {
                        int cz1 = kcz+52-shift_z;
                        for(kcy=0; kcy<ncy; kcy++)
                        {
                            int cy1 = kcy+68-shift_y;
                            if(nboundary(cy1,cz1))
                            {
                                int cy2min = pTable[cz1*ncy_k+cy1].cymin + shift_y;
                                int cz2min = pTable[cz1*ncy_k+cy1].czmin + shift_z;
                                int c2total = pTable[cz1*ncy_k+cy1+1].bias;
                                //if(c2total>0) printf("c2total=%d\n",c2total);
                                // if(pTable[cz1*ncy_k/4+cy1].bias>0) printf("ktotal=%d\n",pTable[cz1*ncy_k/4+cy1].bias);
                                int ktotal;
                                for(ktotal=pTable[cz1*ncy_k+cy1].bias; ktotal<c2total; ktotal++)
                                {//printf("123\n");
                                    int cy2 = pIncre[ktotal]._y+cy2min;
                                    int cz2 = pIncre[ktotal]._z+cz2min;
                                    // printf("cy2= %d cz2= %d ",cy2,cz2);
                                    
                                    if(boundary(cy2-68,cz2-52))
                                    {
                                        pFwd[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)] += pSparse[ktotal] * phantom[fx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j] / pNorm_LOR[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)];  //(float)LOR_initial;
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

void diff(float *pFW, float *pHist, float *pDiff)
{
    int i;
    for(i = 0; i < ncrow * nccol; i++)
    {
        if(pHist[i] >= 0)
        {
            pDiff[i] = abs(pFW[i]-pHist[i]);
            pDiff[i] = 2*pDiff[i];
        }
    }
}

void norm(float *sol, float *norm2, int ite)
{
    int i;
    for(i=0;i<nimgx*nimgy*nimgz;i++)
        norm2[ite] = norm2[ite] + sol[i]*sol[i];
}


void cost_function(float * pFW, float *pHist, float* norm2, float *costfunction, int ite, float gamma)
{
    int i;
    for(i = 0; i < ncrow * nccol; i++)
        costfunction[ite] = costfunction[ite] + pow((pFW[i]-pHist[i]),2);

    printf(" Fidelity term [%d] = %f\n", ite, costfunction[ite]);
    printf("       TV term [%d] = %f\n", ite, gamma*sqrt(norm2[ite]));

    costfunction[ite] = costfunction[ite] + gamma *sqrt(norm2[ite]);
    printf(" cost function [%d] = %lf\n", ite, costfunction[ite]);
}

void backward(float* pDiff, float* pRecon, float* pNorm_voxel)
{
    int rx_plane = x_plane +60;
    int jlayer1 = 0;
    int jlayer2 = 0;
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
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_table_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inTAB = fopen(name,"rb");
            fread(pTable,sizeof(TABLE_ZIP),nccol_k+1,inTAB);
            INCREMENT_ZIP* pIncre = (INCREMENT_ZIP*) calloc(pTable[nccol_k].bias,sizeof(INCREMENT_ZIP));
            float* pSparse = (float*) calloc(pTable[nccol_k].bias,sizeof(float));
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_incre_plane_%d_%s_", x_plane, kernel_name[i*4+j]);
            inINCRE = fopen(name,"rb");
            fread(pIncre,sizeof(INCREMENT_ZIP),pTable[nccol_k].bias,inINCRE);
            
            sprintf(name, "../EM_recon/system matrix/data/sparse_matrix_plane_%d_%s_",x_plane, kernel_name[i*4+j]);
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
                                        pRecon[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j]+=pDiff[(kcz)*ncy*ncrow+(kcy)*ncrow+(cz2-52)*ncy+(cy2-68)]*pSparse[ktotal]/pNorm_voxel[rx_plane*nimgz*nimgy+(iz*4+i)*nimgy+jy*4+j];
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



void proxtv(float* x_n,/*float* u_n,*/ float* r, float *s, float* r2, float* s2,float* dx, float* dy,float* pold, float* qold, float* p, float* q, float* temp_s, float* temp_r, float* sum_rs,float gamma, float* Y_n)
{
    float weights;
    
    float weight_x=1, weight_y=1, mt=1, t, told=1;
    int rx_plane = x_plane +60;
    int i,k;
/////////////////////proxTV iteration//////////////////////////////
    
    // printf("%d",iteration); 
    for (i=0;i<nimgy*nimgz;i++)
    {
        if (i%nimgy==0)
        	s2[i]=s[i];

        else if (i%nimgy==nimgy-1)
        	s2[i]=-s[i-1];

        else
            s2[i]=s[i]-s[i-1];

        s2[i]=s2[i]*weight_x;
    }
    
    for (i=0;i<nimgz;i++)
    {
        for (k=0;k<nimgy;k++)
        {
            if (i==0)
                r2[i*nimgy+k]=r[i*nimgy+k];

            else if (i==nimgz-1)
                r2[i*nimgy+k]=-1*r[(i-1)*(nimgy)+k];

            else
                r2[i*nimgy+k]=r[i*nimgy+k]-r[(i-1)*nimgy+k];
            
            r2[i]=r2[i]*weight_y;
        }
    }
    /* current solution */
    for (i=0;i<nimgy*nimgz;i++)
    {
        x_n[rx_plane*nimgy*nimgz+i] = Y_n[rx_plane*nimgy*nimgz+i] - gamma*(r2[i]+s2[i]);
        //sol[rx_plane*nimgy*nimgz+i]=Y_n[i]-gamma*(r2[i]+s2[i]);
    }
////////////////////////////////////////////////////////////////////
    //update divergence vectors and project
    /* start gradient_op */
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
            max=sum_rs[i];

        weights=max;

        if (max<1)
            weights=1;
    }
    
    for (i=0;i<nimgy*nimgz;i++)
    {
        p[i]=r[i]/weights;
        q[i]=s[i]/weights;
    }
    
    /* FISTA update */
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





int main()
{
    config();

    // pHist_1:nonTOF
    float* pHist_1 = (float*) calloc(nccol*ncrow,sizeof(float));   
    // Activity Reconstruction
    float* pNorm_LOR_1 = (float*) calloc(nccol*ncrow,sizeof(float));
    float* pNorm_voxel_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pFwd_1 = (float*) calloc(nccol*ncrow,sizeof(float));
    float* pBwd_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pDiff_1 = (float*) calloc(nccol*ncrow,sizeof(float));
    float* pRecon_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    
    // FISTA 1 (Activity)
    float * Y_n_1 = (float*)calloc(nimgx*nimgy*nimgz,sizeof(float));    //applying 2D to 3D
    float * r_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * s_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * pold_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * qold_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    // Prox_tv 
    float * r2_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * s2_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * dx_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * dy_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * p_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * q_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_r_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_s_1 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * sum_rs_1 = (float*)calloc(nimgy*nimgz,sizeof(float));


    float* z_n_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* sol_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* x_n_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pFwd_1n = (float*) calloc(nccol*ncrow,sizeof(float));
    float* norm2_1 = (float*) calloc(100,sizeof(float));
    float* costfunction_1 = (float*) calloc(100,sizeof(float));
    // float* phantom_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    // float* u_n_1 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
 
    // pHist_2:TOF
    float* pHist_2 = (float*) calloc(nccol*ncrow,sizeof(float));
    // Attenuation Reconstruction
    float* pNorm_LOR_2 = (float*) calloc(nccol*ncrow,sizeof(float));
    float* pNorm_voxel_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pFwd_2 = (float*) calloc(nccol*ncrow,sizeof(float));
    float* pBwd_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pDiff_2 = (float*) calloc(nccol*ncrow,sizeof(float));
    float* pRecon_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    
    // FISTA 2 (Attenuation)
    float * Y_n_2 = (float*)calloc(nimgx*nimgy*nimgz,sizeof(float));    //applying 2D to 3D
    float * r_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * s_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * pold_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * qold_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    // Prox_tv
    float * r2_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * s2_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * dx_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * dy_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * p_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * q_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_r_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * temp_s_2 = (float*)calloc(nimgy*nimgz,sizeof(float));
    float * sum_rs_2 = (float*)calloc(nimgy*nimgz,sizeof(float));

    float* z_n_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* sol_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* x_n_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    float* pFwd_2n = (float*) calloc(nccol*ncrow,sizeof(float));
    float* costfunction_2 = (float*) calloc(100,sizeof(float));
    float* norm2_2 = (float*) calloc(100,sizeof(float));
    // float* phantom_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));
    // float* u_n_2 = (float*) calloc(nimgx*nimgy*nimgz,sizeof(float));

    // Define parameters
    float tn = 1, tn1;
    int i,j,k;
    float gamma = 0.01;
    int iteration = 20;
    int rx_plane = x_plane + 60;
    int a1 = -675;  int a2 = -4.95;
    int ix,iy,iz,kk,ss;
    
    memset(pFwd_1n, 0, ncrow*nccol * sizeof(float));
    memset(costfunction_1, 0, 100*sizeof(float));
    memset(norm2_1, 0, 100*sizeof(float));
    memset(pFwd_2n, 0, ncrow*nccol * sizeof(float));
    memset(costfunction_2, 0, 100*sizeof(float));
    memset(norm2_2, 0, 100*sizeof(float));
    // printf("stop 01\n");

    printf(" Read data from histogram\n");
    // Read nonTOF data from histogram
    FILE *iHist_1;
    iHist_1 = fopen("../histogram/histogram_nonTOF_new.dat","rb");
    fread(pHist_1,sizeof(float),nccol*ncrow,iHist_1);
    fclose(iHist_1);
    
    // Read TOF data from histogram
    FILE *iHist_2;
    iHist_2 = fopen("../histogram/histogram_TOF_new.dat","rb");
    fread(pHist_2,sizeof(float),nccol*ncrow,iHist_2);
    fclose(iHist_2);
    // printf("stop 02\n");
    
    // Set arrays to zero
    for(ix=0;ix<120;ix++)
    {
        for(i=0;i<nimgy*nimgz;i++)
        {
            // nonTOF (Activity)
            pBwd_1[ix*nimgz*nimgy+i] = 0.0;
            x_n_1[ix*nimgz*nimgy+i] = pHist_1[ix*nimgz*nimgy+i];
            sol_1[ix*nimgz*nimgy+i] = 0.0; 
            // u_n_1[ix*nimgz*nimgy+i]=0.0;

            // TOF (Attenuation)
            pBwd_2[ix*nimgz*nimgy+i] = 0.0;
            x_n_2[ix*nimgz*nimgy+i] = pHist_2[ix*nimgz*nimgy+i];
            sol_2[ix*nimgz*nimgy+i] = 0.0;
            // u_n_2[ix*nimgz*nimgy+i]=0.0;
        }
    }
    
    
    
    ///////////////////////////////// start ADMM ///////////////////////////////////
    int ite = 0, iter;
    float T[3];
    char out_name[300];
    T[0] = s_time();

    /* Initial memory */
    memset(pNorm_LOR_1,0,nccol*ncrow*sizeof(float));
    memset(pNorm_voxel_1,0,nimgx*nimgy*nimgz*sizeof(float));
    memset(pFwd_1,0,nccol*ncrow*sizeof(float));
    memset(pRecon_1,0,nimgx*nimgy*nimgz*sizeof(float));
    
    memset(pNorm_LOR_2,0,nccol*ncrow*sizeof(float));
    memset(pNorm_voxel_2,0,nimgx*nimgy*nimgz*sizeof(float));
    memset(pFwd_2,0,nccol*ncrow*sizeof(float));
    memset(pRecon_2,0,nimgx*nimgy*nimgz*sizeof(float));
    
    /* Stopping Criterions */
    const int MAX_ITER = 100;  //maximum number of iterations
    const float tol = 1e-4;    //relative error tolerance
    /*
    ///////////////////// Activity & Attenuation Reconstruction ///////////////////
    for(x_plane=0;x_plane<3;x_plane++)
    {
        x_plane = x_plane*(-20);
        printf(" norm_LOR_1 : %d\n", x_plane);
        norm_LOR(pNorm_LOR_1);

        printf(" norm_voxel_1 : %d\n", x_plane);
        norm_voxel(pNorm_voxel_1);

        printf(" norm_LOR_2 : %d\n", x_plane);
        norm_LOR(pNorm_LOR_2);

        printf(" norm_voxel_2 : %d\n", x_plane);
        norm_voxel(pNorm_voxel_2);
        x_plane = x_plane/(-20);
    }
    
    
    // Save pNorm_LOR data temperary
    printf("Save data to pNorm_LOR_1\n");
    sprintf(out_name,"pNorm_LOR_1.DAT");
    FILE *iNorm_LOR_1;
    if((iNorm_LOR_1=fopen(out_name,"wb"))==NULL)
    {
        printf("Can not open file\n");
        exit(0);
    }
    fwrite(pNorm_LOR_1,sizeof(float)*nccol*ncrow,1,iNorm_LOR_1);
    fclose(iNorm_LOR_1);

    sprintf(out_name,"pNorm_LOR_2.DAT");
    FILE *iNorm_LOR_2;
    if((iNorm_LOR_2=fopen(out_name,"wb"))==NULL)
    {
        printf("Can not open file\n");
        exit(0);
    }
    fwrite(pNorm_LOR_2,sizeof(float)*nccol*ncrow,1,iNorm_LOR_2);
    fclose(iNorm_LOR_2);
    */
    // Read data from pNorm_LOR
    printf(" Read data from pNorm_LOR\n");
    FILE *iNorm_LOR_1;
    iNorm_LOR_1 = fopen("pNorm_LOR_1.DAT","rb");
    fread(pNorm_LOR_1,sizeof(float),nccol*ncrow,iNorm_LOR_1);
    fclose(iNorm_LOR_1);

    FILE *iNorm_LOR_2;
    iNorm_LOR_2 = fopen("pNorm_LOR_2.DAT","rb");
    fread(pNorm_LOR_2,sizeof(float),nccol*ncrow,iNorm_LOR_2);
    fclose(iNorm_LOR_2);

    /*
    // pBwd = phantom
    for(x_plane=0;x_plane<3;x_plane++)
    {
        x_plane = x_plane*(-20);
        printf(" start forward_1 : %d\n",x_plane);
        forward(pFwd_1,pBwd_1,pNorm_LOR_1);

        printf(" start forward_2 : %d\n",x_plane);
        forward(pFwd_2,pBwd_2,pNorm_LOR_2);
        x_plane = x_plane/(-20);
    }

    
    printf(" start Diff_1\n");
    diff(pFwd_1,pHist_1,pDiff_1);

    printf(" start Diff_2\n");
    diff(pFwd_2,pHist_2,pDiff_2);
    
    for(x_plane=0;x_plane<3;x_plane++)
    {
        x_plane = x_plane*(-20);
        printf(" start backward_1 : %d\n",x_plane);
        backward(pDiff_1,pRecon_1,pNorm_voxel_1);

        printf(" start backward_2 : %d\n",x_plane);
        backward(pDiff_2,pRecon_2,pNorm_voxel_2);
        x_plane = x_plane/(-20);
    }
    */

    /*
    // Save FBP recon data temperary
    sprintf(out_name,"pRecon_1.DAT");
    FILE *pRecon1;
    if((pRecon1=fopen(out_name,"wb"))==NULL)
    {
        printf("Can not open file\n");
        exit(0);
    }
    fwrite(pRecon_1,sizeof(float)*nimgx*nimgy*nimgz,1,pRecon1);
    fclose(pRecon1);

    sprintf(out_name,"pRecon_2.DAT");
    FILE *pRecon2;
    if((pRecon2=fopen(out_name,"wb"))==NULL)
    {
        printf("Can not open file\n");
        exit(0);
    }
    fwrite(pRecon_2,sizeof(float)*nimgx*nimgy*nimgz,1,pRecon2);
    fclose(pRecon2);
    */
    
    // Read data from pRecon
    printf(" Read data from pRecon\n");
    FILE *iRecon_1;
    iRecon_1 = fopen("pRecon_1.DAT","rb");
    fread(pRecon_1,sizeof(float),nimgx*nimgy*nimgz,iRecon_1);
    fclose(iRecon_1);

    // printf("Read data from pRecon\n");
    FILE *iRecon_2;
    iRecon_2 = fopen("pRecon_2.DAT","rb");
    fread(pRecon_2,sizeof(float),nimgx*nimgy*nimgz,iRecon_2);
    fclose(iRecon_2);
    


    T[1] = s_time()-T[0];
    printf("time : %g sec \n",T[1]);
 
   /* Main Loop */
    printf("\n= Main ADMM loop =\n");
    // while (ite < MAX_ITER)
    // for (iter=0;iter<iteration;iter++)
    // {
        
        // for (ite=0;ite<iteration;ite++)    // only run one iteration to calculate cost_function
        // {
        //     printf("\niteration = %d\n",ite+1);

            /* u_n-gamma*f2.grad(u_n) */
            for(x_plane=0;x_plane<3;x_plane++)
            {
                x_plane = x_plane*(-20);
                printf("\n  x_plane = %d\n", x_plane);

                ////////////////////////// FISTA 1 (Activity) ////////////////////////////////
                // printf(" start L1_norm_prox_1 : %d\n",x_plane);
                for (i=0;i<nimgy*nimgz;i++)
                {
                    rx_plane = x_plane + 60;
                    /* applying 2D to 3D Y_n */
                    Y_n_1[rx_plane*nimgy*nimgz+i] = pBwd_1[rx_plane*nimgy*nimgz+i] - gamma*pRecon_1[rx_plane*nimgy*nimgz+i]; 
                    
                    /*initial*/
                    r_1[i] = 0;
                    s_1[i] = 0;
                    pold_1[i] = r_1[i];
                    qold_1[i] = s_1[i];
                }
                
                ///////////////////////// FISTA 2 (Attenuation) /////////////////////////////
                // printf(" start L1_norm_prox_2 : %d\n",x_plane);
                for (i=0;i<nimgy*nimgz;i++)
                {
                    rx_plane = x_plane + 60;
                    /* applying 2D to 3D Y_n */
                    Y_n_2[rx_plane*nimgy*nimgz+i] = pBwd_2[rx_plane*nimgy*nimgz+i] - gamma*pRecon_2[rx_plane*nimgy*nimgz+i];

                    /*initial*/
                    r_2[i] = 0;
                    s_2[i] = 0;
                    pold_2[i] = r_2[i];
                    qold_2[i] = s_2[i];
                }

                for (ite=0;ite<iteration;ite++)    // only run one iteration to calculate cost_function
                {
                    printf("\niteration = %d\n",ite+1);


                /* Activity */
                // printf(" Y_n_2 - z_n_1 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    Y_n_2[i] = Y_n_2[rx_plane*nimgy*nimgz+i] - z_n_1[rx_plane*nimgy*nimgz+i];

                // printf(" x_n_1 = prox_tv(Y_n_2 - z_n_1) : %d\n",x_plane);
                proxtv(x_n_1,/*u_n_1,*/r_1,s_1,r2_1,s2_1,dx_1,dy_1,pold_1,qold_1,p_1,q_1,temp_s_1,temp_r_1,sum_rs_1,gamma,Y_n_2);
                
                // printf(" x_n_1 + z_n_1 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    x_n_1[i] = x_n_1[rx_plane*nimgy*nimgz+i] + z_n_1[rx_plane*nimgy*nimgz+i];
                
                // printf(" Y_n_1 = prox_tv(x_n_1 + z_n_1) : %d\n",x_plane);
                proxtv(Y_n_1,/*u_n_1,*/r_1,s_1,r2_1,s2_1,dx_1,dy_1,pold_1,qold_1,p_1,q_1,temp_s_1,temp_r_1,sum_rs_1,gamma,x_n_1);
                
                // printf(" update z_n_1 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    z_n_1[i] = z_n_1[rx_plane*nimgy*nimgz+i] - Y_n_1[rx_plane*nimgy*nimgz+i] + x_n_1[rx_plane*nimgy*nimgz+i];

                // printf(" update sol_1: %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    sol_1[i] = a1*x_n_1[rx_plane*nimgy*nimgz+i];
                

                /* Attenuation */
                // printf(" Y_n_1 - z_n_2 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    Y_n_1[i] = Y_n_1[rx_plane*nimgy*nimgz+i] - z_n_2[rx_plane*nimgy*nimgz+i];

                // printf(" x_n_2 = prox_tv(Y_n_1 - z_n_2) : %d\n",x_plane);
                proxtv(x_n_2,/*u_n_2,*/r_2,s_2,r2_2,s2_2,dx_2,dy_2,pold_2,qold_2,p_2,q_2,temp_s_2,temp_r_2,sum_rs_2,gamma,Y_n_1);
                
                // printf(" x_n_2 + z_n_2 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    x_n_2[i] = x_n_2[rx_plane*nimgy*nimgz+i] + z_n_2[rx_plane*nimgy*nimgz+i];
                
                // printf(" Y_n_2 = prox_tv(x_n_2 + z_n_2) : %d\n",x_plane);
                proxtv(Y_n_2,/*u_n_1,*/r_2,s_2,r2_2,s2_2,dx_2,dy_2,pold_2,qold_2,p_2,q_2,temp_s_2,temp_r_2,sum_rs_2,gamma,x_n_2);
                
                // printf(" update z_n_2 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    z_n_2[i] = z_n_2[rx_plane*nimgy*nimgz+i] - Y_n_2[rx_plane*nimgy*nimgz+i] + x_n_2[rx_plane*nimgy*nimgz+i];

                // printf(" update sol_2 : %d\n",x_plane);
                for(i=0;i<nimgy*nimgz;i++)
                    sol_2[i] = a2*x_n_2[rx_plane*nimgy*nimgz+i];


                // if(ite==iteration)
                // {
                    /* Activity */
                    if(ite<20)
                    // if(((ite+1)%5==0) || (ite<5))   //only check iterations: 1~5, 10, 15, ...
                    {
                        if(x_plane==-40)
                        {
                            forward(pFwd_1n,sol_1,pNorm_LOR_1);
                            norm(sol_1,norm2_1,ite);
                            
                            printf(" norm2_1 = %f\n",norm2_1[ite]);
                            cost_function(pFwd_1n,pHist_1,norm2_1,costfunction_1,ite,gamma);
                            // cost_function(sol_1,pHist_1,norm2_1,costfunction_1,ite,gamma);
                            
                            sprintf(out_name,"DHAPET_ite%d_activity.DAT",ite+1);
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


                    /* Attenuation */
                    if(ite<20)
                    // if(((ite+1)%5==0) || (ite<5))
                    {
                        if(x_plane==-40)
                        {
                            forward(pFwd_2n,sol_2,pNorm_LOR_2);
                            norm(sol_2,norm2_2,ite);

                            printf("\n norm2_2 = %f\n",norm2_2[ite]);
                            cost_function(pFwd_2n,pHist_2,norm2_2,costfunction_2,ite,gamma);
                            // cost_function(sol_2,pHist_1,norm2_1,costfunction_1,ite,gamma);
                            
                            sprintf(out_name,"DHAPET_ite%d_attenuation.DAT",ite+1);
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
                }

                x_plane = x_plane/(-20);
            // }

            // if (costfunction_1[ite] < tol)
            //     break;

            // if (costfunction_2[ite] < tol)
            //     break;
                // ite++;
            
        }
        
    // }


    T[2] = s_time()-T[1];
    printf("time: %g sec \n",T[2]);
    
    
/////////////////////////////////stop ADMM/////////////////////////////////////////////
    
    /* Clear memory */
    free(norm2_1);
    free(x_n_1);
    // free(u_n_1);
    free(sol_1);
    // free(phantom_1);
    free(pFwd_1);
    free(pBwd_1);
    free(pHist_1);
    free(pDiff_1);
    free(pRecon_1);
    free(pNorm_LOR_1);
    free(pNorm_voxel_1);
    free(costfunction_1);
    
    free(norm2_2);
    free(x_n_2);
    // free(u_n_2);
    free(sol_2);
    // free(phantom_2);
    free(pFwd_2);
    free(pBwd_2);
    free(pHist_2);
    free(pDiff_2);
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
